#include "smesh_sideset.hpp"
#include "smesh_adjacency.hpp"
#include "smesh_alloc.hpp"
#include "smesh_file_extensions.hpp"
#include "smesh_glob.hpp"
#include "smesh_mesh.hpp"
#include "smesh_read.hpp"
#include "smesh_semistructured.hpp"
#include "smesh_sidesets.hpp"
#include "smesh_sshex8_graph.hpp"
#include "smesh_sspyramid_graph.hpp"
#include "smesh_ssquad4_graph.hpp"
#include "smesh_ssquad4_mesh.hpp"
#include "smesh_sstet4.hpp"
#include "smesh_sstet4_graph.hpp"
#include "smesh_sswedge_graph.hpp"
#include "smesh_tracer.hpp"
#include "smesh_write.hpp"

#ifdef SMESH_ENABLE_MPI
#include "smesh_alltoallv.impl.hpp"
#include "smesh_decompose.hpp"
#include "smesh_distributed_base.hpp"
#include "smesh_distributed_write.hpp"
#endif

#include <algorithm>
#include <cstddef>
#include <cstring>
#include <fstream>
#include <iostream>
#include <list>
#include <unordered_map>
#include <vector>

namespace smesh {

    class Sideset::Impl final {
    public:
        std::shared_ptr<Communicator>          comm;
        std::shared_ptr<Buffer<element_idx_t>> parent;
        std::shared_ptr<Buffer<i16>>           lfi;
        block_idx_t                            block_id{0};
        SharedBuffer<large_idx_t>              element_mapping;
    };

    Sideset::Sideset(const std::shared_ptr<Communicator>          &comm,
                     const std::shared_ptr<Buffer<element_idx_t>> &parent,
                     const std::shared_ptr<Buffer<i16>>           &lfi,
                     block_idx_t                                   block_id,
                     const SharedBuffer<large_idx_t>              &element_mapping)
        : impl_(std::make_unique<Impl>()) {
        impl_->comm            = comm;
        impl_->parent          = parent;
        impl_->lfi             = lfi;
        impl_->block_id        = block_id;
        impl_->element_mapping = element_mapping;
    }

    Sideset::Sideset() : impl_(std::make_unique<Impl>()) {}
    Sideset::~Sideset() = default;

    ptrdiff_t Sideset::size() const { return impl_->parent->size(); }

    std::shared_ptr<Communicator> Sideset::comm() const { return impl_->comm; }

    std::shared_ptr<Sideset> Sideset::create(const std::shared_ptr<Communicator>          &comm,
                                             const std::shared_ptr<Buffer<element_idx_t>> &parent,
                                             const std::shared_ptr<Buffer<i16>>           &lfi,
                                             block_idx_t                                   block_id,
                                             const SharedBuffer<large_idx_t>              &element_mapping) {
        return std::make_shared<Sideset>(comm, parent, lfi, block_id, element_mapping);
    }

    std::shared_ptr<Sideset> Sideset::create_from_file(const std::shared_ptr<Communicator> &comm,
                                                       const Path                          &path,
                                                       block_idx_t                          block_id) {
        auto ret = std::make_shared<Sideset>();
        if (ret->read(comm, path, block_id) != SMESH_SUCCESS) return nullptr;
        return ret;
    }

    int Sideset::write(const Path &path) const {
        if (!comm()->rank()) {
            smesh::create_directory(path.to_string());
        }

        Path parent_path = path / ("parent." + std::string(TypeToString<element_idx_t>::value()));
        Path lfi_path    = path / ("lfi." + std::string(TypeToString<i16>::value()));

#ifdef SMESH_ENABLE_MPI
        comm()->barrier();
        if (comm()->size() > 1) {
            u64 n_sides_local  = impl_->parent->size();
            u64 n_sides_global = n_sides_local;
            MPI_Allreduce(MPI_IN_PLACE, &n_sides_global, 1, mpi_type<u64>(), MPI_SUM, comm()->get());

            if (impl_->element_mapping) {
                auto parent   = create_host_buffer<element_idx_t>(n_sides_local);
                auto b_parent = parent->data();
                for (u64 i = 0; i < n_sides_local; i++) {
                    b_parent[i] = impl_->element_mapping->data()[impl_->parent->data()[i]];
                }
                array_write_convert_from_extension(comm()->get(), parent_path, b_parent, n_sides_local, n_sides_global);
            } else {
                array_write_convert_from_extension(
                        comm()->get(), parent_path, impl_->parent->data(), impl_->parent->size(), n_sides_global);
            }
            array_write_convert_from_extension(comm()->get(), lfi_path, impl_->lfi->data(), impl_->lfi->size(), n_sides_global);
        } else
#endif  // SMESH_ENABLE_MPI
        {
            array_write(parent_path, impl_->parent->data(), impl_->parent->size());
            array_write(lfi_path, impl_->lfi->data(), impl_->lfi->size());
        }

        if (!comm()->rank()) {
            std::ofstream os((path / "meta.yaml").to_string());

            if (os.good()) {
                os << "# Automatically generated by smesh_Sideset.cpp\n";
                os << "size: " << impl_->parent->size() << "\n";
                os << "parent: parent.raw\n";
                os << "lfi: lfi.int16.raw\n";
                os << "rpath: true\n";
            } else {
                SMESH_ERROR("Unable to open sideset meta.yaml file\n");
                return SMESH_FAILURE;
            }

            os.close();
        }

        return SMESH_SUCCESS;
    }

    static std::vector<std::shared_ptr<Sideset>> filter_sidesets_to_owned(
            const std::shared_ptr<Mesh>                          &mesh,
            std::vector<std::shared_ptr<Sideset>>                 sidesets) {
        if (mesh->comm()->size() == 1) {
            return sidesets;
        }

        for (auto &sideset : sidesets) {
            auto block = mesh->block(sideset->block_id());
            const ptrdiff_t n_owned = block->n_elements_owned();

            const ptrdiff_t n_sides = sideset->size();
            auto            parent  = sideset->parent()->data();
            auto            lfi     = sideset->lfi()->data();

            std::vector<element_idx_t> owned_parent;
            std::vector<i16>           owned_lfi;
            owned_parent.reserve(static_cast<size_t>(n_sides));
            owned_lfi.reserve(static_cast<size_t>(n_sides));

            for (ptrdiff_t i = 0; i < n_sides; ++i) {
                if (parent[i] < n_owned) {
                    owned_parent.push_back(parent[i]);
                    owned_lfi.push_back(lfi[i]);
                }
            }

            element_idx_t *new_parent = nullptr;
            i16           *new_lfi    = nullptr;
            if (!owned_parent.empty()) {
                new_parent = (element_idx_t *)SMESH_ALLOC(owned_parent.size() * sizeof(element_idx_t));
                new_lfi    = (i16 *)SMESH_ALLOC(owned_lfi.size() * sizeof(i16));
                std::memcpy(new_parent, owned_parent.data(), owned_parent.size() * sizeof(element_idx_t));
                std::memcpy(new_lfi, owned_lfi.data(), owned_lfi.size() * sizeof(i16));
            }

            sideset = std::make_shared<Sideset>(mesh->comm(),
                                                smesh::manage_host_buffer(
                                                    static_cast<ptrdiff_t>(owned_parent.size()), new_parent),
                                                smesh::manage_host_buffer(
                                                    static_cast<ptrdiff_t>(owned_lfi.size()), new_lfi),
                                                sideset->block_id(),
                                                block->element_mapping());
        }

        return sidesets;
    }

    std::vector<std::shared_ptr<Sideset>> Sideset::create_from_selector(
            const std::shared_ptr<Mesh>                                         &mesh,
            const std::function<bool(const geom_t, const geom_t, const geom_t)> &selector,
            const std::vector<std::string>                                      &block_names) {
        SMESH_TRACE_SCOPE("Sideset::create_from_selector");

        const int dim    = mesh->spatial_dimension();
        auto      points = mesh->points()->data();

        size_t                                n_blocks = mesh->n_blocks();
        std::vector<std::shared_ptr<Sideset>> sidesets;

        for (size_t b = 0; b < n_blocks; b++) {
            auto block = mesh->block(b);
            if (!block_names.empty() &&  //
                std::find(block_names.begin(), block_names.end(), block->name()) == block_names.end()) {
                continue;
            }

            const enum ElemType element_type = block->element_type();
            const bool          is_ss        = is_semistructured_type(element_type);
            const enum ElemType family       = is_ss ? ss_source_family(element_type) : element_type;
            int                 corners[8];
            int                 n_corners = 0;
            if (is_ss && !ss_source_family_corners(family, semistructured_level(element_type), corners, &n_corners)) {
                SMESH_ERROR("Sideset::create_from_selector: SS family %s is not supported\n", type_to_string(family));
                return {};
            }
            (void)n_corners;

            const int      ns = elem_num_sides(family);
            LocalSideTable lst;
            lst.fill(family);

            const ptrdiff_t nelements = block->n_elements();
            auto            elements  = block->elements()->data();

#if 0
    std::list<element_idx_t> parent_list;
    std::list<i16> lfi_list;
    for (ptrdiff_t e = 0; e < nelements; e++) {
      for (int s = 0; s < ns; s++) {
        // Barycenter of face
        double p[3] = {0, 0, 0};

        for (int ln = 0; ln < nnxs; ln++) {
          const idx_t node = elements[lst(s, ln)][e];

          for (int d = 0; d < dim; d++) {
            p[d] += points[d][node];
          }
        }

        for (int d = 0; d < dim; d++) {
          p[d] /= nnxs;
        }

        if (selector(p[0], p[1], p[2])) {
          parent_list.push_back(e);
          lfi_list.push_back(s);
        }
      }
    }

    const ptrdiff_t nparents = parent_list.size();
    auto parent = create_host_buffer<element_idx_t>(nparents);
    auto lfi = create_host_buffer<i16>(nparents);

    {
      ptrdiff_t idx = 0;
      for (auto p : parent_list) {
        parent->data()[idx++] = p;
      }
    }

    {
      ptrdiff_t idx = 0;
      for (auto l : lfi_list) {
        lfi->data()[idx++] = l;
      }
    }

#else
            auto      selected   = create_host_buffer<u8>(nelements * ns);
            auto      b_selected = selected->data();
            ptrdiff_t nparents   = 0;
#pragma omp parallel for reduction(+ : nparents)
            for (ptrdiff_t e = 0; e < nelements; e++) {
                for (int s = 0; s < ns; s++) {
                    // Barycenter of face
                    double p[3] = {0, 0, 0};

                    for (int ln = 0; ln < lst.nnxs_side[s]; ln++) {
                        const int   soa_row = is_ss ? corners[lst(s, ln)] : lst(s, ln);
                        const idx_t node    = elements[soa_row][e];

                        for (int d = 0; d < dim; d++) {
                            p[d] += points[d][node];
                        }
                    }

                    for (int d = 0; d < dim; d++) {
                        p[d] /= lst.nnxs_side[s];
                    }

                    if (selector(p[0], p[1], p[2])) {
                        b_selected[e * ns + s] = 1;
                        nparents++;
                    }
                }
            }

            auto parent = create_host_buffer<element_idx_t>(nparents);
            auto lfi    = create_host_buffer<i16>(nparents);

            auto b_parent = parent->data();
            auto b_lfi    = lfi->data();

            ptrdiff_t offset = 0;
            // #pragma omp parallel for
            for (ptrdiff_t e = 0; e < nelements; e++) {
                for (int s = 0; s < ns; s++) {
                    if (b_selected[e * ns + s]) {
                        b_parent[offset] = e;
                        b_lfi[offset]    = s;
                        offset++;
                    }
                }
            }

#endif

            sidesets.push_back(std::make_shared<Sideset>(
                    mesh->comm(),
                    parent,
                    lfi,
                    static_cast<block_idx_t>(b),
                    mesh->comm()->size() > 1 ? block->element_mapping() : nullptr));
        }

        return filter_sidesets_to_owned(mesh, sidesets);
    }

    std::vector<std::shared_ptr<Sideset>> Sideset::create_from_batch_selector(
            const std::shared_ptr<Mesh> &mesh,
            const std::function<
                    void(const ptrdiff_t, const geom_t *const, const geom_t *const, const geom_t *const, u8 *const selected)>
                                           &selector,
            const std::vector<std::string> &block_names) {
        SMESH_TRACE_SCOPE("Sideset::create_from_batch_selector");

        const int dim    = mesh->spatial_dimension();
        auto      points = mesh->points()->data();

        size_t                                n_blocks = mesh->n_blocks();
        std::vector<std::shared_ptr<Sideset>> sidesets;

        for (size_t b = 0; b < n_blocks; b++) {
            auto block = mesh->block(b);
            if (!block_names.empty() &&  //
                std::find(block_names.begin(), block_names.end(), block->name()) == block_names.end()) {
                continue;
            }

            const enum ElemType element_type = block->element_type();
            const bool          is_ss        = is_semistructured_type(element_type);
            const enum ElemType family       = is_ss ? ss_source_family(element_type) : element_type;
            int                 corners[8];
            int                 n_corners = 0;
            if (is_ss && !ss_source_family_corners(family, semistructured_level(element_type), corners, &n_corners)) {
                SMESH_ERROR("Sideset::create_from_batch_selector: SS family %s is not supported\n", type_to_string(family));
                return {};
            }
            (void)n_corners;

            const int      ns = elem_num_sides(family);
            LocalSideTable lst;
            lst.fill(family);

            const ptrdiff_t nelements = block->n_elements();
            auto            elements  = block->elements()->data();

            auto      selected   = create_host_buffer<u8>(nelements * ns);
            auto      b_selected = selected->data();
            ptrdiff_t nparents   = 0;

            const ptrdiff_t nelements_per_batch = 16;
            const ptrdiff_t ns_per_batch        = nelements_per_batch * ns;

#pragma omp parallel
            {
                auto fx        = create_host_buffer<geom_t>(ns_per_batch);
                auto fy        = create_host_buffer<geom_t>(ns_per_batch);
                auto fz        = create_host_buffer<geom_t>(ns_per_batch);
                auto fselected = create_host_buffer<u8>(ns_per_batch);

                auto b_fx        = fx->data();
                auto b_fy        = fy->data();
                auto b_fz        = fz->data();
                auto b_fselected = fselected->data();

#pragma omp for schedule(static) reduction(+ : nparents)
                for (ptrdiff_t e = 0; e < nelements; e += nelements_per_batch) {
                    const ptrdiff_t batch_size = std::min(nelements_per_batch, (nelements - e));

                    for (ptrdiff_t b = 0; b < batch_size; b++) {
                        for (int s = 0; s < ns; s++) {
                            geom_t p[3] = {0, 0, 0};

                            for (int ln = 0; ln < lst.nnxs_side[s]; ln++) {
                                const int   soa_row = is_ss ? corners[lst(s, ln)] : lst(s, ln);
                                const idx_t node    = elements[soa_row][e + b];

                                for (int d = 0; d < dim; d++) {
                                    p[d] += points[d][node];
                                }
                            }

                            for (int d = 0; d < dim; d++) {
                                p[d] /= lst.nnxs_side[s];
                            }

                            const ptrdiff_t idx = b * ns + s;
                            b_fx[idx]           = p[0];
                            b_fy[idx]           = p[1];
                            b_fz[idx]           = p[2];
                            b_fselected[idx]    = 0;
                        }
                    }

                    selector(batch_size * ns, b_fx, b_fy, b_fz, b_fselected);

                    for (ptrdiff_t b = 0; b < batch_size; b++) {
                        for (int s = 0; s < ns; s++) {
                            const ptrdiff_t idx = b * ns + s;
                            if (b_fselected[idx]) {
                                b_selected[(e + b) * ns + s] = 1;
                                nparents++;
                            }
                        }
                    }
                }
            }

            auto parent = create_host_buffer<element_idx_t>(nparents);
            auto lfi    = create_host_buffer<i16>(nparents);

            auto b_parent = parent->data();
            auto b_lfi    = lfi->data();

            ptrdiff_t offset = 0;
            for (ptrdiff_t e = 0; e < nelements; e++) {
                for (int s = 0; s < ns; s++) {
                    if (b_selected[e * ns + s]) {
                        b_parent[offset] = e;
                        b_lfi[offset]    = s;
                        offset++;
                    }
                }
            }

            sidesets.push_back(std::make_shared<Sideset>(
                    mesh->comm(),
                    parent,
                    lfi,
                    static_cast<block_idx_t>(b),
                    mesh->comm()->size() > 1 ? block->element_mapping() : nullptr));
        }

        return filter_sidesets_to_owned(mesh, sidesets);
    }

    std::vector<std::shared_ptr<Sideset>> Sideset::create_from_plane(const std::shared_ptr<Mesh>    &mesh,
                                                                     const geom_t                    normal_x,
                                                                     const geom_t                    normal_y,
                                                                     const geom_t                    normal_z,
                                                                     const geom_t                    distance,
                                                                     const geom_t                    tol,
                                                                     const std::vector<std::string> &block_names) {
        SMESH_TRACE_SCOPE("Sideset::create_from_plane");
        return create_from_batch_selector(
                mesh,
                [=](const ptrdiff_t                    batch_size,
                    const geom_t *const SMESH_RESTRICT x,
                    const geom_t *const SMESH_RESTRICT y,
                    const geom_t *const SMESH_RESTRICT z,
                    u8 *const SMESH_RESTRICT           selected) {
#pragma omp simd
                    for (ptrdiff_t i = 0; i < batch_size; i++) {
                        selected[i] = std::abs(normal_x * x[i] + normal_y * y[i] + normal_z * z[i] - distance) < tol;
                    }
                },
                block_names);
    }

    int Sideset::read(const std::shared_ptr<Communicator> &comm, const Path &folder, block_idx_t block_id) {
        SMESH_TRACE_SCOPE("Sideset::read");
        if (comm->size() > 1) {
            SMESH_ERROR("Sideset::read is not supported for distributed runs\n");
            return SMESH_FAILURE;
        }

        impl_->comm = comm;
        ptrdiff_t      nlocal{0}, ncheck{0};
        element_idx_t *parent{nullptr};
        i16           *lfi{nullptr};

        auto parent_file = detect_files(folder / "parent.*", {"raw", "int16", "int32", "int64"});
        auto lfi_file    = detect_files(folder / "lfi.*", {"raw", "int16", "int32", "int64"});

        if (parent_file.empty() || lfi_file.empty()) {
            SMESH_ERROR("Unable to find parent or lfi file in sideset at %s\n", folder.c_str());
            return SMESH_FAILURE;
        }

        int ret = SMESH_SUCCESS;
        ret |= array_read_convert_from_extension(parent_file[0], &parent, &nlocal);
        ret |= array_read_convert_from_extension(lfi_file[0], &lfi, &ncheck);

        impl_->parent = smesh::manage_host_buffer(nlocal, parent);
        impl_->lfi    = smesh::manage_host_buffer(nlocal, lfi);

        if (ncheck != nlocal) {
            SMESH_ERROR("Inconsistent array sizes in sideset at %s\n", folder.c_str());
            ret = SMESH_FAILURE;
        }

        impl_->block_id = block_id;
        return ret;
    }

    int Sideset::read_and_redistibute(const std::shared_ptr<Mesh> &mesh, const Path &path, block_idx_t block_id) {
#ifdef SMESH_ENABLE_MPI
        const int rank = mesh->comm()->rank();
        const int size = mesh->comm()->size();
        if (size > 1) {
            i64 n_sides = 0;
            if (rank == 0) {
                if (read(Communicator::self(), path, block_id) != SMESH_SUCCESS) {
                    n_sides = -1;
                } else {
                    n_sides = static_cast<i64>(this->size());
                }
            }
            mesh->comm()->broadcast(&n_sides, 1, 0);
            if (n_sides < 0) {
                return SMESH_FAILURE;
            }

            mesh->comm()->broadcast(&block_id, 1, 0);
            impl_->block_id = block_id;
            impl_->comm     = mesh->comm();

            element_idx_t *parent = nullptr;
            i16           *lfi    = nullptr;
            if (n_sides > 0) {
                parent = (element_idx_t *)SMESH_ALLOC(static_cast<size_t>(n_sides) * sizeof(element_idx_t));
                lfi    = (i16 *)SMESH_ALLOC(static_cast<size_t>(n_sides) * sizeof(i16));
                if (rank == 0) {
                    std::memcpy(parent, impl_->parent->data(), static_cast<size_t>(n_sides) * sizeof(element_idx_t));
                    std::memcpy(lfi, impl_->lfi->data(), static_cast<size_t>(n_sides) * sizeof(i16));
                }
                mesh->comm()->broadcast(parent, static_cast<int>(n_sides), 0);
                mesh->comm()->broadcast(lfi, static_cast<int>(n_sides), 0);
            }

            impl_->parent = smesh::manage_host_buffer(static_cast<ptrdiff_t>(n_sides), parent);
            impl_->lfi    = smesh::manage_host_buffer(static_cast<ptrdiff_t>(n_sides), lfi);
            return redistribute(mesh);
        }
#endif
        if (read(mesh->comm(), path, block_id) != SMESH_SUCCESS) {
            return SMESH_FAILURE;
        }
        return redistribute(mesh);
    }

    int Sideset::redistribute(const std::shared_ptr<Mesh> &mesh) {
        if (mesh->comm()->size() == 1) {
            return SMESH_SUCCESS;
        }

#ifdef SMESH_ENABLE_MPI
        auto block = mesh->block(impl_->block_id);
        if (!block) {
            SMESH_ERROR("Sideset::redistribute: invalid block_id %d\n", impl_->block_id);
            return SMESH_FAILURE;
        }

        const ptrdiff_t n_owned = block->n_elements_owned();
        auto            elem_map = block->element_mapping();
        if (n_owned > 0 && !elem_map) {
            SMESH_ERROR("Sideset::redistribute: block %s has no element_mapping\n", block->name().c_str());
            return SMESH_FAILURE;
        }

        large_idx_t n_global = 0;
        if (n_owned > 0) {
            for (ptrdiff_t i = 0; i < n_owned; ++i) {
                n_global = std::max(n_global, elem_map->data()[i] + 1);
            }
        }
        n_global = mesh->comm()->max(n_global);

        std::vector<element_idx_t> invert(static_cast<size_t>(n_global),
                                          invalid_idx<element_idx_t>());
        for (ptrdiff_t i = 0; i < n_owned; ++i) {
            const large_idx_t serial_id = elem_map->data()[i];
            if (serial_id >= 0 && serial_id < n_global) {
                invert[static_cast<size_t>(serial_id)] = static_cast<element_idx_t>(i);
            }
        }

        const ptrdiff_t n_sides = size();
        auto            parent  = impl_->parent->data();
        auto            lfi     = impl_->lfi->data();

        std::vector<element_idx_t> local_parent;
        std::vector<i16>           local_lfi;
        local_parent.reserve(static_cast<size_t>(n_sides));
        local_lfi.reserve(static_cast<size_t>(n_sides));

        for (ptrdiff_t i = 0; i < n_sides; ++i) {
            const large_idx_t serial_parent = static_cast<large_idx_t>(parent[i]);
            if (serial_parent >= 0 && serial_parent < n_global) {
                const element_idx_t local_e = invert[static_cast<size_t>(serial_parent)];
                if (local_e != invalid_idx<element_idx_t>()) {
                    local_parent.push_back(local_e);
                    local_lfi.push_back(lfi[i]);
                }
            }
        }

        element_idx_t *new_parent =
            (element_idx_t *)SMESH_ALLOC(local_parent.size() * sizeof(element_idx_t));
        i16 *new_lfi = (i16 *)SMESH_ALLOC(local_lfi.size() * sizeof(i16));
        if (!local_parent.empty()) {
            std::memcpy(new_parent, local_parent.data(), local_parent.size() * sizeof(element_idx_t));
            std::memcpy(new_lfi, local_lfi.data(), local_lfi.size() * sizeof(i16));
        }

        impl_->parent          = smesh::manage_host_buffer(static_cast<ptrdiff_t>(local_parent.size()), new_parent);
        impl_->lfi             = smesh::manage_host_buffer(static_cast<ptrdiff_t>(local_lfi.size()), new_lfi);
        impl_->element_mapping = block->element_mapping();

        return SMESH_SUCCESS;
#endif

        SMESH_ERROR("Sideset::redistribute is not supported for distributed runs\n");
        return SMESH_FAILURE;
    }

    std::shared_ptr<Buffer<element_idx_t>> Sideset::parent() const { return impl_->parent; }
    std::shared_ptr<Buffer<i16>>           Sideset::lfi() const { return impl_->lfi; }
    block_idx_t                            Sideset::block_id() const { return impl_->block_id; }

    std::shared_ptr<Buffer<idx_t>> create_nodeset_from_sidesets(const std::shared_ptr<Mesh>                 &mesh,
                                                                const std::vector<std::shared_ptr<Sideset>> &sidesets) {
        if (sidesets.empty()) {
            return nullptr;
        }

        if (sidesets.size() == 1) {
            return create_nodeset_from_sideset(mesh, sidesets[0]);
        }

        bool any_ss = false;
        for (const auto &ss : sidesets) {
            if (is_semistructured_type(mesh->element_type(ss->block_id()))) {
                any_ss = true;
                break;
            }
        }

        // extract_nodeset_from_sidesets uses LocalSideTable, which only knows HEX8/HEX27/
        // PROTEUS_HEX8 corners — not PROTEUS_HEX27/64/125/... Full SS faces go through
        // sshex8_extract_nodeset_from_sideset (create_nodeset_from_sideset).
        if (any_ss) {
            std::vector<idx_t> ids;
            for (const auto &ss : sidesets) {
                auto ns = create_nodeset_from_sideset(mesh, ss);
                if (!ns || ns->size() == 0) {
                    continue;
                }
                auto d = ns->data();
                ids.insert(ids.end(), d, d + ns->size());
            }
            std::sort(ids.begin(), ids.end());
            ids.erase(std::unique(ids.begin(), ids.end()), ids.end());

            const ptrdiff_t n     = static_cast<ptrdiff_t>(ids.size());
            idx_t          *nodes = nullptr;
            if (n > 0) {
                nodes = (idx_t *)SMESH_ALLOC(static_cast<size_t>(n) * sizeof(idx_t));
                std::memcpy(nodes, ids.data(), static_cast<size_t>(n) * sizeof(idx_t));
            }
            return smesh::manage_host_buffer(n, nodes);
        }

        block_idx_t                  n_sidesets = sidesets.size();
        std::vector<enum ElemType>   element_type(n_sidesets);
        std::vector<idx_t **>        elems(n_sidesets);
        std::vector<ptrdiff_t>       n_surf_elements(n_sidesets);
        std::vector<element_idx_t *> parent_element(n_sidesets);
        std::vector<i16 *>           side_idx(n_sidesets);

        for (block_idx_t k = 0; k < n_sidesets; k++) {
            auto ss    = sidesets[k];
            auto block = mesh->block(ss->block_id());

            element_type[k]    = block->element_type();
            elems[k]           = block->elements()->data();
            n_surf_elements[k] = ss->parent()->size();
            parent_element[k]  = ss->parent()->data();
            side_idx[k]        = ss->lfi()->data();
        }

        ptrdiff_t n_nodes_out{0};
        idx_t    *nodes_out{nullptr};

        extract_nodeset_from_sidesets(n_sidesets,
                                      element_type.data(),
                                      elems.data(),
                                      n_surf_elements.data(),
                                      parent_element.data(),
                                      side_idx.data(),
                                      &n_nodes_out,
                                      &nodes_out);

        return smesh::manage_host_buffer(n_nodes_out, nodes_out);
    }

    // std::shared_ptr<Buffer<idx_t>>
    // create_nodeset_from_sideset(const std::shared_ptr<Mesh> &mesh,
    //                             const std::shared_ptr<Sideset> &sideset) {

    // ptrdiff_t n_nodes{0};
    // idx_t *nodes{nullptr};
    //   if (extract_nodeset_from_sideset(
    //           mesh->element_type(), mesh->elements()->data(),
    //           sideset->parent()->size(), sideset->parent()->data(),
    //           sideset->lfi()->data(), &n_nodes, &nodes) != SMESH_SUCCESS) {
    //     SMESH_ERROR("Unable to extract nodeset from sideset!\n");
    //   }

    //   return smesh::manage_host_buffer(n_nodes, nodes);
    // }

    static enum ElemType sstet4_surface_type(const int L) {
        switch (L) {
            case 1:
                return TRISHELL3;
            case 2:
                return TRISHELL6;
            case 3:
                return TRISHELL10;
            default:
                return INVALID;
        }
    }

    std::pair<enum ElemType, std::shared_ptr<Buffer<idx_t *>>> create_surface_from_sideset_semistructured(
            const std::shared_ptr<Mesh>    &ssmesh,
            const std::shared_ptr<Sideset> &sideset) {
        auto block = ssmesh->block(sideset->block_id());
        if (!block) {
            SMESH_ERROR("Unable to find block %d for sideset surface extract!\n", static_cast<int>(sideset->block_id()));
            return {INVALID, nullptr};
        }

        const enum ElemType et = block->element_type();
        if (!is_semistructured_type(et)) {
            SMESH_ERROR("Mesh is not a semistructured mesh!\n");
            return {INVALID, nullptr};
        }

        const enum ElemType family = ss_source_family(et);
        const int           L      = semistructured_level(et);
        const ptrdiff_t     n_surf = sideset->parent()->size();

        if (family == HEX8) {
            auto ss_sides = smesh::create_host_buffer<idx_t>((L + 1) * (L + 1), n_surf);
            if (n_surf > 0 && sshex8_extract_surface_from_sideset(L,
                                                                  block->elements()->data(),
                                                                  n_surf,
                                                                  sideset->parent()->data(),
                                                                  sideset->lfi()->data(),
                                                                  ss_sides->data()) != SMESH_SUCCESS) {
                SMESH_ERROR("Unable to extract surface from sideset!\n");
            }
            return {shell_type(side_type(et)), ss_sides};
        }

        if (family == TET4) {
            auto ss_sides = smesh::create_host_buffer<idx_t>(sstet4_n_tri(L), n_surf);
            if (n_surf > 0 && sstet4_extract_surface_from_sideset(L,
                                                                  block->elements()->data(),
                                                                  n_surf,
                                                                  sideset->parent()->data(),
                                                                  sideset->lfi()->data(),
                                                                  ss_sides->data()) != SMESH_SUCCESS) {
                SMESH_ERROR("Unable to extract surface from sideset!\n");
            }
            return {sstet4_surface_type(L), ss_sides};
        }

        if (family == QUAD4) {
            auto ss_sides = smesh::create_host_buffer<idx_t>(L + 1, n_surf);
            if (n_surf > 0 && ssquad4_extract_surface_from_sideset(L,
                                                                   block->elements()->data(),
                                                                   n_surf,
                                                                   sideset->parent()->data(),
                                                                   sideset->lfi()->data(),
                                                                   ss_sides->data()) != SMESH_SUCCESS) {
                SMESH_ERROR("Unable to extract surface from sideset!\n");
            }
            return {shell_type(side_type(et)), ss_sides};
        }

        if (family == WEDGE6 || family == PYRAMID5) {
            if (n_surf <= 0) {
                return {INVALID, nullptr};
            }
            const int  first        = sideset->lfi()->data()[0];
            const bool is_quad_side = (family == WEDGE6) ? (first < 3) : (first >= 4);
            for (ptrdiff_t i = 0; i < n_surf; ++i) {
                const int s = sideset->lfi()->data()[i];
                const bool q = (family == WEDGE6) ? (s < 3) : (s >= 4);
                if (q != is_quad_side) {
                    SMESH_ERROR("create_surface_from_sideset_semistructured: mixed tri/quad sideset arity\n");
                    return {INVALID, nullptr};
                }
            }
            if (is_quad_side) {
                auto ss_sides = smesh::create_host_buffer<idx_t>((L + 1) * (L + 1), n_surf);
                int  err      = SMESH_SUCCESS;
                if (family == WEDGE6) {
                    err = sswedge_extract_surface_from_sideset(L,
                                                               block->elements()->data(),
                                                               n_surf,
                                                               sideset->parent()->data(),
                                                               sideset->lfi()->data(),
                                                               ss_sides->data());
                } else {
                    err = sspyramid_extract_surface_from_sideset(L,
                                                                 block->elements()->data(),
                                                                 n_surf,
                                                                 sideset->parent()->data(),
                                                                 sideset->lfi()->data(),
                                                                 ss_sides->data());
                }
                if (err != SMESH_SUCCESS) {
                    SMESH_ERROR("Unable to extract surface from sideset!\n");
                }
                return {shell_type(proteus_quad_type(L)), ss_sides};
            }
            auto ss_sides = smesh::create_host_buffer<idx_t>(sstet4_n_tri(L), n_surf);
            int  err      = SMESH_SUCCESS;
            if (family == WEDGE6) {
                err = sswedge_extract_surface_from_sideset(L,
                                                           block->elements()->data(),
                                                           n_surf,
                                                           sideset->parent()->data(),
                                                           sideset->lfi()->data(),
                                                           ss_sides->data());
            } else {
                err = sspyramid_extract_surface_from_sideset(L,
                                                             block->elements()->data(),
                                                             n_surf,
                                                             sideset->parent()->data(),
                                                             sideset->lfi()->data(),
                                                             ss_sides->data());
            }
            if (err != SMESH_SUCCESS) {
                SMESH_ERROR("Unable to extract surface from sideset!\n");
            }
            return {sstet4_surface_type(L), ss_sides};
        }

        SMESH_ERROR("create_surface_from_sideset_semistructured: family %s is not implemented\n", type_to_string(family));
        return {INVALID, nullptr};
    }

    std::pair<enum ElemType, std::shared_ptr<Buffer<idx_t *>>> create_surface_from_sideset(
            const std::shared_ptr<Mesh>    &mesh,
            const std::shared_ptr<Sideset> &sideset) {
        if (is_semistructured_type(mesh->element_type(sideset->block_id()))) {
            return create_surface_from_sideset_semistructured(mesh, sideset);
        }

        auto block = mesh->block(sideset->block_id());
        const ptrdiff_t n_surf = sideset->parent()->size();
        enum ElemType   face_st = INVALID;
        int             nnxs    = 0;
        if (n_surf > 0) {
            face_st = shell_type(side_type(block->element_type(), sideset->lfi()->data()[0]));
            nnxs    = elem_num_nodes(face_st);
        } else if (elem_sides_homogeneous(block->element_type())) {
            face_st = shell_type(side_type(block->element_type()));
            nnxs    = elem_num_nodes(face_st);
        } else {
            return {INVALID, nullptr};
        }

        auto surface = smesh::create_host_buffer<idx_t>(nnxs, n_surf);

        if (extract_surface_from_sideset(block->element_type(),
                                         block->elements()->data(),
                                         sideset->parent()->size(),
                                         sideset->parent()->data(),
                                         sideset->lfi()->data(),
                                         surface->data()) != SMESH_SUCCESS) {
            SMESH_ERROR("Unable to create surface from sideset!");
        }

        // printf("Type: %s->%s\n", type_to_string(block->element_type()),
        // type_to_string(st)); surface->print();

        return {face_st, surface};
    }

    std::shared_ptr<Buffer<idx_t>> create_nodeset_from_sideset_semistructured(const std::shared_ptr<Mesh>    &ss,
                                                                              const std::shared_ptr<Sideset> &sideset) {
        auto block = ss->block(sideset->block_id());
        if (!block) {
            SMESH_ERROR("Unable to find block %d for sideset nodeset extract!\n", static_cast<int>(sideset->block_id()));
            return nullptr;
        }

        const enum ElemType et = block->element_type();
        if (!is_semistructured_type(et)) {
            SMESH_ERROR("Mesh is not a semistructured mesh!\n");
            return nullptr;
        }

        const enum ElemType family = ss_source_family(et);
        const int           L      = semistructured_level(et);
        const ptrdiff_t     n_surf = sideset->parent()->size();
        ptrdiff_t           n_nodes{0};
        idx_t              *nodes{nullptr};

        if (n_surf == 0) {
            return smesh::manage_host_buffer(n_nodes, nodes);
        }

        if (family == HEX8) {
            SMESH_TRACE_SCOPE("sshex8_extract_nodeset_from_sideset");
            if (sshex8_extract_nodeset_from_sideset(L,
                                                    block->elements()->data(),
                                                    sideset->parent()->size(),
                                                    sideset->parent()->data(),
                                                    sideset->lfi()->data(),
                                                    &n_nodes,
                                                    &nodes) != SMESH_SUCCESS) {
                SMESH_ERROR("Unable to extract nodeset from sideset!\n");
            }
            return smesh::manage_host_buffer(n_nodes, nodes);
        }

        if (family == TET4) {
            SMESH_TRACE_SCOPE("sstet4_extract_nodeset_from_sideset");
            if (sstet4_extract_nodeset_from_sideset(L,
                                                    block->elements()->data(),
                                                    sideset->parent()->size(),
                                                    sideset->parent()->data(),
                                                    sideset->lfi()->data(),
                                                    &n_nodes,
                                                    &nodes) != SMESH_SUCCESS) {
                SMESH_ERROR("Unable to extract nodeset from sideset!\n");
            }
            return smesh::manage_host_buffer(n_nodes, nodes);
        }

        if (family == QUAD4) {
            SMESH_TRACE_SCOPE("ssquad4_extract_nodeset_from_sideset");
            if (ssquad4_extract_nodeset_from_sideset(L,
                                                     block->elements()->data(),
                                                     sideset->parent()->size(),
                                                     sideset->parent()->data(),
                                                     sideset->lfi()->data(),
                                                     &n_nodes,
                                                     &nodes) != SMESH_SUCCESS) {
                SMESH_ERROR("Unable to extract nodeset from sideset!\n");
            }
            return smesh::manage_host_buffer(n_nodes, nodes);
        }

        if (family == WEDGE6) {
            SMESH_TRACE_SCOPE("sswedge_extract_nodeset_from_sideset");
            if (sswedge_extract_nodeset_from_sideset(L,
                                                     block->elements()->data(),
                                                     sideset->parent()->size(),
                                                     sideset->parent()->data(),
                                                     sideset->lfi()->data(),
                                                     &n_nodes,
                                                     &nodes) != SMESH_SUCCESS) {
                SMESH_ERROR("Unable to extract nodeset from sideset!\n");
            }
            return smesh::manage_host_buffer(n_nodes, nodes);
        }

        if (family == PYRAMID5) {
            SMESH_TRACE_SCOPE("sspyramid_extract_nodeset_from_sideset");
            if (sspyramid_extract_nodeset_from_sideset(L,
                                                       block->elements()->data(),
                                                       sideset->parent()->size(),
                                                       sideset->parent()->data(),
                                                       sideset->lfi()->data(),
                                                       &n_nodes,
                                                       &nodes) != SMESH_SUCCESS) {
                SMESH_ERROR("Unable to extract nodeset from sideset!\n");
            }
            return smesh::manage_host_buffer(n_nodes, nodes);
        }

        SMESH_ERROR("create_nodeset_from_sideset_semistructured: family %s is not implemented\n", type_to_string(family));
        return nullptr;
    }

    // Element-owned sidesets miss constrained nodes whose owner rank only
    // holds the incident element as a ghost. Broadcast GIDs so every rank
    // that stores the node (owned or halo) includes it in the nodeset.
    static std::shared_ptr<Buffer<idx_t>> synchronize_nodeset_gids(const std::shared_ptr<Mesh>          &mesh,
                                                                   const std::shared_ptr<Buffer<idx_t>> &local) {
#ifdef SMESH_ENABLE_MPI
        if (!mesh || !mesh->comm() || mesh->comm()->size() <= 1 || !mesh->is_distributed() || !mesh->distributed()) {
            return local;
        }
        auto nmap_buf = mesh->distributed()->node_mapping();
        if (!nmap_buf || nmap_buf->size() == 0) {
            return local;
        }

        const large_idx_t *const nmap           = nmap_buf->data();
        const ptrdiff_t          n_local_nodes  = mesh->n_nodes();
        std::vector<large_idx_t> send;
        if (local && local->size() > 0) {
            send.reserve(static_cast<size_t>(local->size()));
            auto d = local->data();
            for (ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(local->size()); ++i) {
                const idx_t n = d[i];
                if (n >= 0 && static_cast<ptrdiff_t>(n) < n_local_nodes) {
                    send.push_back(nmap[n]);
                }
            }
            std::sort(send.begin(), send.end());
            send.erase(std::unique(send.begin(), send.end()), send.end());
        }

        MPI_Comm comm  = mesh->comm()->get();
        const int nsend  = static_cast<int>(send.size());
        const int nranks = mesh->comm()->size();
        std::vector<int> counts(static_cast<size_t>(nranks), 0);
        std::vector<int> displs(static_cast<size_t>(nranks), 0);
        SMESH_MPI_CATCH(MPI_Allgather(&nsend, 1, MPI_INT, counts.data(), 1, MPI_INT, comm));
        int ntotal = 0;
        for (int r = 0; r < nranks; ++r) {
            displs[static_cast<size_t>(r)] = ntotal;
            ntotal += counts[static_cast<size_t>(r)];
        }
        std::vector<large_idx_t> recv(static_cast<size_t>(ntotal));
        SMESH_MPI_CATCH(MPI_Allgatherv(send.empty() ? nullptr : send.data(),
                                       nsend,
                                       mpi_type<large_idx_t>(),
                                       recv.empty() ? nullptr : recv.data(),
                                       counts.data(),
                                       displs.data(),
                                       mpi_type<large_idx_t>(),
                                       comm));

        std::unordered_map<large_idx_t, idx_t> gid_to_local;
        gid_to_local.reserve(static_cast<size_t>(n_local_nodes));
        for (ptrdiff_t i = 0; i < n_local_nodes; ++i) {
            gid_to_local.emplace(nmap[i], static_cast<idx_t>(i));
        }

        std::vector<idx_t> out;
        out.reserve(local ? static_cast<size_t>(local->size()) : 0);
        for (int i = 0; i < ntotal; ++i) {
            auto it = gid_to_local.find(recv[static_cast<size_t>(i)]);
            if (it != gid_to_local.end()) {
                out.push_back(it->second);
            }
        }
        std::sort(out.begin(), out.end());
        out.erase(std::unique(out.begin(), out.end()), out.end());

        const ptrdiff_t n     = static_cast<ptrdiff_t>(out.size());
        idx_t          *nodes = nullptr;
        if (n > 0) {
            nodes = (idx_t *)SMESH_ALLOC(static_cast<size_t>(n) * sizeof(idx_t));
            std::memcpy(nodes, out.data(), static_cast<size_t>(n) * sizeof(idx_t));
        }
        return manage_host_buffer(n, nodes);
#else
        (void)mesh;
        return local;
#endif
    }

    std::shared_ptr<Buffer<idx_t>> create_nodeset_from_sideset(const std::shared_ptr<Mesh>    &mesh,
                                                               const std::shared_ptr<Sideset> &sideset) {
        std::shared_ptr<Buffer<idx_t>> local;
        if (is_semistructured_type(mesh->element_type(sideset->block_id()))) {
            local = create_nodeset_from_sideset_semistructured(mesh, sideset);
        } else {
            auto block = mesh->block(sideset->block_id());

            ptrdiff_t n_nodes{0};
            idx_t    *nodes{nullptr};
            if (extract_nodeset_from_sideset(block->element_type(),
                                             block->elements()->data(),
                                             sideset->parent()->size(),
                                             sideset->parent()->data(),
                                             sideset->lfi()->data(),
                                             &n_nodes,
                                             &nodes) != SMESH_SUCCESS) {
                SMESH_ERROR("Unable to extract nodeset from sideset!\n");
            }

            local = smesh::manage_host_buffer(n_nodes, nodes);
        }
        return synchronize_nodeset_gids(mesh, local);
    }

    std::pair<enum ElemType, std::shared_ptr<Buffer<idx_t *>>> create_surface_from_sidesets(
            const std::shared_ptr<Mesh>                 &mesh,
            const std::vector<std::shared_ptr<Sideset>> &sidesets) {
        if (sidesets.empty()) {
            return {INVALID, nullptr};
        }

        if (sidesets.size() != 1) {
            SMESH_ERROR(
                    "create_surface_from_sidesets is not supported for multiple "
                    "sidesets!\n");
            return {INVALID, nullptr};
        }

        if (is_semistructured_type(mesh->element_type(sidesets[0]->block_id()))) {
            return create_surface_from_sideset_semistructured(mesh, sidesets[0]);
        } else {
            return create_surface_from_sideset(mesh, sidesets[0]);
        }
    }

    void Sideset::print(std::ostream &os) const {
        os << "Sideset: " << block_id() << "\n";
        os << "Size: " << parent()->size() << "\n";

        parent()->print(os);
        lfi()->print(os);
    }

}  // namespace smesh


