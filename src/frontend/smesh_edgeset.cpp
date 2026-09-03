#include "smesh_edgeset.hpp"
#include "smesh_adjacency.hpp"
#include "smesh_alloc.hpp"
#include "smesh_edgesets.hpp"
#include "smesh_elem_type.hpp"
#include "smesh_file_extensions.hpp"
#include "smesh_glob.hpp"
#include "smesh_mesh.hpp"
#include "smesh_read.hpp"
#include "smesh_tracer.hpp"
#include "smesh_write.hpp"

#ifdef SMESH_ENABLE_MPI
#include "smesh_distributed_base.hpp"
#include "smesh_distributed_write.hpp"
#endif

#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace smesh {

    class Edgeset::Impl final {
    public:
        std::shared_ptr<Communicator>          comm;
        std::shared_ptr<Buffer<element_idx_t>> parent;
        std::shared_ptr<Buffer<i16>>           lei;
        block_idx_t                            block_id{0};
        SharedBuffer<large_idx_t>              element_mapping;
    };

    Edgeset::Edgeset(const std::shared_ptr<Communicator>          &comm,
                     const std::shared_ptr<Buffer<element_idx_t>> &parent,
                     const std::shared_ptr<Buffer<i16>>           &lei,
                     block_idx_t                                   block_id,
                     const SharedBuffer<large_idx_t>              &element_mapping)
        : impl_(std::make_unique<Impl>()) {
        impl_->comm            = comm;
        impl_->parent          = parent;
        impl_->lei             = lei;
        impl_->block_id        = block_id;
        impl_->element_mapping = element_mapping;
    }

    Edgeset::Edgeset() : impl_(std::make_unique<Impl>()) {}
    Edgeset::~Edgeset() = default;

    ptrdiff_t Edgeset::size() const { return impl_->parent ? (ptrdiff_t)impl_->parent->size() : 0; }

    std::shared_ptr<Communicator> Edgeset::comm() const { return impl_->comm; }

    std::shared_ptr<Edgeset> Edgeset::create(const std::shared_ptr<Communicator>          &comm,
                                             const std::shared_ptr<Buffer<element_idx_t>> &parent,
                                             const std::shared_ptr<Buffer<i16>>           &lei,
                                             block_idx_t                                   block_id,
                                             const SharedBuffer<large_idx_t>              &element_mapping) {
        return std::make_shared<Edgeset>(comm, parent, lei, block_id, element_mapping);
    }

    std::shared_ptr<Edgeset> Edgeset::create_from_file(const std::shared_ptr<Communicator> &comm,
                                                       const Path                          &path,
                                                       block_idx_t                          block_id) {
        auto ret = std::make_shared<Edgeset>();
        if (ret->read(comm, path, block_id) != SMESH_SUCCESS) {
            return nullptr;
        }
        return ret;
    }

    static int parse_set_meta(const Path &folder,
                              bool       *has_size,
                              ptrdiff_t  *meta_size,
                              bool       *has_block_id,
                              block_idx_t *meta_block_id) {
        *has_size     = false;
        *has_block_id = false;
        const Path meta = folder / "meta.yaml";
        if (!meta.exists()) {
            return SMESH_SUCCESS;
        }

        std::ifstream ifs(meta.to_string());
        if (!ifs.good()) {
            SMESH_ERROR("Unable to open edgeset meta.yaml file\n");
            return SMESH_FAILURE;
        }

        std::string line;
        while (std::getline(ifs, line)) {
            const auto hash = line.find('#');
            if (hash != std::string::npos) {
                line.resize(hash);
            }
            const auto start = line.find_first_not_of(" \t");
            if (start == std::string::npos) {
                continue;
            }
            line = line.substr(start);
            if (line.compare(0, 5, "size:") == 0) {
                *has_size  = true;
                *meta_size = static_cast<ptrdiff_t>(std::strtoll(line.c_str() + 5, nullptr, 10));
            } else if (line.compare(0, 9, "block_id:") == 0) {
                *has_block_id  = true;
                *meta_block_id = static_cast<block_idx_t>(std::strtol(line.c_str() + 9, nullptr, 10));
            }
        }
        return SMESH_SUCCESS;
    }

    int Edgeset::write(const Path &path) const {
        if (!comm()->rank()) {
            smesh::create_directory(path.to_string());
        }

        Path parent_path = path / ("parent." + std::string(TypeToString<element_idx_t>::value()));
        Path lei_path    = path / ("lei." + std::string(TypeToString<i16>::value()));

        u64 n_edges_global = static_cast<u64>(size());

#ifdef SMESH_ENABLE_MPI
        comm()->barrier();
        if (comm()->size() > 1) {
            const u64 n_edges_local = n_edges_global;
            MPI_Allreduce(MPI_IN_PLACE, &n_edges_global, 1, mpi_type<u64>(), MPI_SUM, comm()->get());

            if (impl_->element_mapping) {
                auto parent   = create_host_buffer<element_idx_t>(n_edges_local);
                auto b_parent = parent->data();
                const large_idx_t   *b_map   = impl_->element_mapping->data();
                const element_idx_t *b_local = impl_->parent->data();
                for (u64 i = 0; i < n_edges_local; i++) {
                    b_parent[i] = b_map[b_local[i]];
                }
                array_write_convert_from_extension(
                        comm()->get(), parent_path, b_parent, n_edges_local, n_edges_global);
            } else {
                array_write_convert_from_extension(
                        comm()->get(), parent_path, impl_->parent->data(), impl_->parent->size(), n_edges_global);
            }
            array_write_convert_from_extension(
                    comm()->get(), lei_path, impl_->lei->data(), impl_->lei->size(), n_edges_global);
        } else
#endif
        {
            array_write(parent_path, impl_->parent->data(), impl_->parent->size());
            array_write(lei_path, impl_->lei->data(), impl_->lei->size());
        }

        if (!comm()->rank()) {
            std::ofstream os((path / "meta.yaml").to_string());
            if (os.good()) {
                os << "# Automatically generated by smesh_edgeset.cpp\n";
                os << "size: " << n_edges_global << "\n";
                os << "block_id: " << static_cast<long long>(impl_->block_id) << "\n";
                os << "parent: parent." << TypeToString<element_idx_t>::value() << "\n";
                os << "lei: lei." << TypeToString<i16>::value() << "\n";
                os << "rpath: true\n";
            } else {
                SMESH_ERROR("Unable to open edgeset meta.yaml file\n");
                return SMESH_FAILURE;
            }
        }
        return SMESH_SUCCESS;
    }

    int Edgeset::read(const std::shared_ptr<Communicator> &comm, const Path &folder, block_idx_t block_id) {
        SMESH_TRACE_SCOPE("Edgeset::read");
        if (comm->size() > 1) {
            SMESH_ERROR("Edgeset::read is not supported for distributed runs\n");
            return SMESH_FAILURE;
        }

        impl_->comm = comm;
        ptrdiff_t      nlocal{0}, ncheck{0};
        element_idx_t *parent{nullptr};
        i16           *lei{nullptr};

        auto parent_file = detect_files(folder / "parent.*", {"raw", "int16", "int32", "int64"});
        auto lei_file    = detect_files(folder / "lei.*", {"raw", "int16", "int32", "int64"});

        if (parent_file.empty() || lei_file.empty()) {
            SMESH_ERROR("Unable to find parent or lei file in edgeset at %s\n", folder.c_str());
            return SMESH_FAILURE;
        }

        int ret = SMESH_SUCCESS;
        ret |= array_read_convert_from_extension(parent_file[0], &parent, &nlocal);
        ret |= array_read_convert_from_extension(lei_file[0], &lei, &ncheck);

        impl_->parent = smesh::manage_host_buffer(nlocal, parent);
        impl_->lei    = smesh::manage_host_buffer(nlocal, lei);

        if (ncheck != nlocal) {
            SMESH_ERROR("Inconsistent array sizes in edgeset at %s\n", folder.c_str());
            ret = SMESH_FAILURE;
        }

        bool        has_size      = false;
        bool        has_block_id  = false;
        ptrdiff_t   meta_size     = 0;
        block_idx_t meta_block_id = 0;
        if (parse_set_meta(folder, &has_size, &meta_size, &has_block_id, &meta_block_id) != SMESH_SUCCESS) {
            return SMESH_FAILURE;
        }
        if (has_size && meta_size != nlocal) {
            SMESH_ERROR("Edgeset meta size %td does not match array length %td at %s\n",
                        meta_size,
                        nlocal,
                        folder.c_str());
            ret = SMESH_FAILURE;
        }

        impl_->block_id = has_block_id ? meta_block_id : block_id;
        return ret;
    }

    int Edgeset::read_and_redistibute(const std::shared_ptr<Mesh> &mesh, const Path &path, block_idx_t block_id) {
#ifdef SMESH_ENABLE_MPI
        const int rank = mesh->comm()->rank();
        const int size = mesh->comm()->size();
        if (size > 1) {
            i64 n_edges = 0;
            if (rank == 0) {
                if (read(Communicator::self(), path, block_id) != SMESH_SUCCESS) {
                    n_edges = -1;
                } else {
                    n_edges  = static_cast<i64>(this->size());
                    block_id = impl_->block_id;
                }
            }
            mesh->comm()->broadcast(&n_edges, 1, 0);
            if (n_edges < 0) {
                return SMESH_FAILURE;
            }

            mesh->comm()->broadcast(&block_id, 1, 0);
            impl_->block_id = block_id;
            impl_->comm     = mesh->comm();

            auto parent_buf = create_host_buffer<element_idx_t>((size_t)n_edges);
            auto lei_buf    = create_host_buffer<i16>((size_t)n_edges);
            if (rank == 0 && n_edges > 0) {
                std::memcpy(parent_buf->data(),
                            impl_->parent->data(),
                            (size_t)n_edges * sizeof(element_idx_t));
                std::memcpy(lei_buf->data(), impl_->lei->data(), (size_t)n_edges * sizeof(i16));
            }
            if (n_edges > 0) {
                mesh->comm()->broadcast(parent_buf->data(), (int)n_edges, 0);
                mesh->comm()->broadcast(lei_buf->data(), (int)n_edges, 0);
            }

            impl_->parent = parent_buf;
            impl_->lei    = lei_buf;
            return redistribute(mesh);
        }
#endif
        if (read(mesh->comm(), path, block_id) != SMESH_SUCCESS) {
            return SMESH_FAILURE;
        }
        return redistribute(mesh);
    }

    int Edgeset::redistribute(const std::shared_ptr<Mesh> &mesh) {
        if (mesh->comm()->size() == 1) {
            return SMESH_SUCCESS;
        }

#ifdef SMESH_ENABLE_MPI
        auto block = mesh->block(impl_->block_id);
        if (!block) {
            SMESH_ERROR("Edgeset::redistribute: invalid block_id %d\n", impl_->block_id);
            return SMESH_FAILURE;
        }

        const ptrdiff_t n_owned  = block->n_elements_owned();
        auto            elem_map = block->element_mapping();
        if (n_owned > 0 && !elem_map) {
            SMESH_ERROR("Edgeset::redistribute: block %s has no element_mapping\n", block->name().c_str());
            return SMESH_FAILURE;
        }

        large_idx_t        n_global = 0;
        const large_idx_t *d_emap   = (n_owned > 0) ? elem_map->data() : nullptr;
        if (n_owned > 0) {
            for (ptrdiff_t i = 0; i < n_owned; ++i) {
                const large_idx_t gid = d_emap[i] + 1;
                if (gid > n_global) {
                    n_global = gid;
                }
            }
        }
        n_global = mesh->comm()->max(n_global);

        if (n_global == 0) {
            impl_->parent          = create_host_buffer<element_idx_t>(0);
            impl_->lei             = create_host_buffer<i16>(0);
            impl_->element_mapping = block->element_mapping();
            return SMESH_SUCCESS;
        }

        auto           invert = create_host_buffer<element_idx_t>((size_t)n_global);
        element_idx_t *d_inv  = invert->data();
        for (large_idx_t i = 0; i < n_global; ++i) {
            d_inv[i] = invalid_idx<element_idx_t>();
        }
        for (ptrdiff_t i = 0; i < n_owned; ++i) {
            const large_idx_t serial_id = d_emap[i];
            if (serial_id >= 0 && serial_id < n_global) {
                d_inv[serial_id] = (element_idx_t)i;
            }
        }

        const ptrdiff_t      n_edges = size();
        const element_idx_t *parent  = impl_->parent->data();
        const i16           *lei     = impl_->lei->data();

        ptrdiff_t n_keep = 0;
        for (ptrdiff_t i = 0; i < n_edges; ++i) {
            const large_idx_t serial_parent = (large_idx_t)parent[i];
            if (serial_parent >= 0 && serial_parent < n_global &&
                d_inv[serial_parent] != invalid_idx<element_idx_t>()) {
                ++n_keep;
            }
        }

        auto           new_parent = create_host_buffer<element_idx_t>((size_t)n_keep);
        auto           new_lei    = create_host_buffer<i16>((size_t)n_keep);
        element_idx_t *d_np       = n_keep > 0 ? new_parent->data() : nullptr;
        i16           *d_nl       = n_keep > 0 ? new_lei->data() : nullptr;
        ptrdiff_t      w          = 0;
        for (ptrdiff_t i = 0; i < n_edges; ++i) {
            const large_idx_t serial_parent = (large_idx_t)parent[i];
            if (serial_parent >= 0 && serial_parent < n_global) {
                const element_idx_t local_e = d_inv[serial_parent];
                if (local_e != invalid_idx<element_idx_t>()) {
                    d_np[w] = local_e;
                    d_nl[w] = lei[i];
                    ++w;
                }
            }
        }

        impl_->parent          = new_parent;
        impl_->lei             = new_lei;
        impl_->element_mapping = block->element_mapping();
        return SMESH_SUCCESS;
#endif

        SMESH_ERROR("Edgeset::redistribute is not supported for distributed runs\n");
        return SMESH_FAILURE;
    }

    std::shared_ptr<Buffer<element_idx_t>> Edgeset::parent() const { return impl_->parent; }
    std::shared_ptr<Buffer<i16>>           Edgeset::lei() const { return impl_->lei; }
    block_idx_t                            Edgeset::block_id() const { return impl_->block_id; }
    SharedBuffer<large_idx_t>              Edgeset::element_mapping() const { return impl_->element_mapping; }

    int Edgeset::remap_parents(const element_idx_t *old_to_new, ptrdiff_t n) {
        if (!old_to_new || n < 0 || !impl_->parent) {
            SMESH_ERROR("Edgeset::remap_parents: invalid arguments\n");
            return SMESH_FAILURE;
        }
        if (impl_->parent->mem_space() != MEMORY_SPACE_HOST) {
            SMESH_ERROR("Edgeset::remap_parents requires host parent buffer\n");
            return SMESH_FAILURE;
        }
        auto            p   = impl_->parent->data();
        const ptrdiff_t nes = size();
        for (ptrdiff_t i = 0; i < nes; ++i) {
            const element_idx_t old = p[i];
            if (old < 0 || static_cast<ptrdiff_t>(old) >= n) {
                SMESH_ERROR("Edgeset::remap_parents: parent %td out of range [0, %td)\n",
                            static_cast<ptrdiff_t>(old),
                            n);
                return SMESH_FAILURE;
            }
            p[i] = old_to_new[old];
        }
        return SMESH_SUCCESS;
    }

    int Edgeset::apply_element_permutation(const element_idx_t *new_to_old, ptrdiff_t n) {
        if (!new_to_old || n < 0) {
            SMESH_ERROR("Edgeset::apply_element_permutation: invalid arguments\n");
            return SMESH_FAILURE;
        }
        if (n == 0) {
            return SMESH_SUCCESS;
        }
        auto           old_to_new = create_host_buffer<element_idx_t>((size_t)n);
        element_idx_t *d_otn      = old_to_new->data();
        for (ptrdiff_t i = 0; i < n; ++i) {
            d_otn[i] = invalid_idx<element_idx_t>();
        }
        for (ptrdiff_t neu = 0; neu < n; ++neu) {
            const element_idx_t old = new_to_old[neu];
            if (old < 0 || (ptrdiff_t)old >= n) {
                SMESH_ERROR("Edgeset::apply_element_permutation: gather index %td out of range\n",
                            (ptrdiff_t)old);
                return SMESH_FAILURE;
            }
            d_otn[old] = (element_idx_t)neu;
        }
        for (ptrdiff_t old = 0; old < n; ++old) {
            if (d_otn[old] == invalid_idx<element_idx_t>()) {
                SMESH_ERROR("Edgeset::apply_element_permutation: permutation is not a bijection\n");
                return SMESH_FAILURE;
            }
        }
        return remap_parents(d_otn, n);
    }

    int remap_edgesets(const std::vector<std::shared_ptr<Edgeset>> &edgesets,
                       block_idx_t                                  block_id,
                       const element_idx_t                         *old_to_new,
                       ptrdiff_t                                    n) {
        for (size_t i = 0; i < edgesets.size(); ++i) {
            const auto &es = edgesets[i];
            if (!es || es->block_id() != block_id) {
                continue;
            }
            if (es->remap_parents(old_to_new, n) != SMESH_SUCCESS) {
                return SMESH_FAILURE;
            }
        }
        return SMESH_SUCCESS;
    }

    std::pair<enum ElemType, std::shared_ptr<Buffer<idx_t *>>> create_edges_from_edgeset(
            const std::shared_ptr<Mesh>    &mesh,
            const std::shared_ptr<Edgeset> &edgeset) {
        if (!mesh || !edgeset) {
            return {INVALID, nullptr};
        }
        auto block = mesh->block(edgeset->block_id());
        if (!block) {
            SMESH_ERROR("create_edges_from_edgeset: invalid block_id %d\n", (int)edgeset->block_id());
            return {INVALID, nullptr};
        }
        const enum ElemType et   = block->element_type();
        LocalEdgeTable      let;
        if (let.fill(et) != SMESH_SUCCESS) {
            LocalEdgeTable::report_unsupported("create_edges_from_edgeset", et);
            return {INVALID, nullptr};
        }
        const int       nne = let.nnxe > 0 ? let.nnxe : 2;
        const ptrdiff_t n   = edgeset->size();
        auto            soa = create_host_buffer<idx_t>(nne, n);
        if (n > 0) {
            if (extract_edges_from_edgeset(et,
                                           block->elements()->data(),
                                           n,
                                           edgeset->parent()->data(),
                                           edgeset->lei()->data(),
                                           soa->data()) != SMESH_SUCCESS) {
                return {INVALID, nullptr};
            }
        }
        const enum ElemType edge_et = (nne == 3) ? EDGE3 : EDGE2;
        return {edge_et, soa};
    }

    void Edgeset::print(std::ostream &os) const {
        os << "Edgeset: " << block_id() << "\n";
        os << "Size: " << size() << "\n";
        if (parent()) {
            parent()->print(os);
        }
        if (lei()) {
            lei()->print(os);
        }
    }

}  // namespace smesh
