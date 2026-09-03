#include "smesh_sideset.hpp"
#include "smesh_adjacency.hpp"
#include "smesh_alloc.hpp"
#include "smesh_file_extensions.hpp"
#include "smesh_glob.hpp"
#include "smesh_mesh.hpp"
#include "smesh_read.hpp"
#include "smesh_semistructured.hpp"
#include "smesh_sidesets.hpp"
#include "smesh_sort.hpp"
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

#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
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

    static int parse_sideset_meta(const Path &folder,
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
            SMESH_ERROR("Unable to open sideset meta.yaml file\n");
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

    int Sideset::write(const Path &path) const {
        if (!comm()->rank()) {
            smesh::create_directory(path.to_string());
        }

        Path parent_path = path / ("parent." + std::string(TypeToString<element_idx_t>::value()));
        Path lfi_path    = path / ("lfi." + std::string(TypeToString<i16>::value()));

        u64 n_sides_global = static_cast<u64>(impl_->parent->size());

#ifdef SMESH_ENABLE_MPI
        comm()->barrier();
        if (comm()->size() > 1) {
            const u64 n_sides_local = n_sides_global;
            MPI_Allreduce(MPI_IN_PLACE, &n_sides_global, 1, mpi_type<u64>(), MPI_SUM, comm()->get());

            if (impl_->element_mapping) {
                auto parent   = create_host_buffer<element_idx_t>(n_sides_local);
                auto b_parent = parent->data();
                const large_idx_t   *b_map    = impl_->element_mapping->data();
                const element_idx_t *b_local  = impl_->parent->data();
                for (u64 i = 0; i < n_sides_local; i++) {
                    b_parent[i] = b_map[b_local[i]];
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
                os << "size: " << n_sides_global << "\n";
                os << "block_id: " << static_cast<long long>(impl_->block_id) << "\n";
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

    static std::shared_ptr<Sideset> sideset_from_buffers(const std::shared_ptr<Mesh>           &mesh,
                                                         const block_idx_t                      block_id,
                                                         const SharedBuffer<element_idx_t>     &parent,
                                                         const SharedBuffer<i16>               &lfi,
                                                         const SharedBuffer<large_idx_t>       &mapping) {
        return std::make_shared<Sideset>(mesh->comm(), parent, lfi, block_id, mapping);
    }

    static int block_name_selected(const std::vector<std::string> &block_names, const std::string &name) {
        if (block_names.empty()) {
            return 1;
        }
        for (size_t i = 0; i < block_names.size(); ++i) {
            if (block_names[i] == name) {
                return 1;
            }
        }
        return 0;
    }

    static std::vector<std::shared_ptr<Sideset>> filter_sidesets_to_owned(
            const std::shared_ptr<Mesh>           &mesh,
            std::vector<std::shared_ptr<Sideset>>  sidesets) {
        if (mesh->comm()->size() == 1) {
            return sidesets;
        }

        for (size_t k = 0; k < sidesets.size(); ++k) {
            auto                &sideset = sidesets[k];
            auto                 block   = mesh->block(sideset->block_id());
            const ptrdiff_t      n_owned = block->n_elements_owned();
            const ptrdiff_t      n_sides = sideset->size();
            const element_idx_t *parent  = sideset->parent()->data();
            const i16           *lfi     = sideset->lfi()->data();

            ptrdiff_t n_keep = 0;
            for (ptrdiff_t i = 0; i < n_sides; ++i) {
                n_keep += (parent[i] < n_owned);
            }
            if (n_keep == n_sides) {
                continue;
            }

            auto            new_parent = create_host_buffer<element_idx_t>((size_t)n_keep);
            auto            new_lfi    = create_host_buffer<i16>((size_t)n_keep);
            element_idx_t  *d_np       = n_keep > 0 ? new_parent->data() : nullptr;
            i16            *d_nl       = n_keep > 0 ? new_lfi->data() : nullptr;
            ptrdiff_t       w          = 0;
            for (ptrdiff_t i = 0; i < n_sides; ++i) {
                if (parent[i] < n_owned) {
                    d_np[w] = parent[i];
                    d_nl[w] = lfi[i];
                    ++w;
                }
            }
            sideset = sideset_from_buffers(mesh, sideset->block_id(), new_parent, new_lfi, block->element_mapping());
        }

        return sidesets;
    }

    static void fill_side_node_key(const bool            is_ss,
                                   const int            *corners,
                                   const LocalSideTable &lst,
                                   idx_t *const *const   elems,
                                   const element_idx_t   e,
                                   const i16             s,
                                   idx_t *const          key,
                                   int                  *nn) {
        *nn = lst.nnxs_side[s];
        for (int ln = 0; ln < *nn; ++ln) {
            const int soa_row = is_ss ? corners[lst(s, ln)] : lst(s, ln);
            key[ln]           = elems[soa_row][e];
        }
        sort(key, (size_t)*nn);
        for (int ln = *nn; ln < LocalSideTable::MAX_NUM_NODES_PER_SIDE; ++ln) {
            key[ln] = invalid_idx<idx_t>();
        }
    }

    /// Keep one copy of each unique face (sorted corner node ids). Owner is
    /// (block_id, parent, lfi) lexicographic — block index ascending, then local
    /// parent, then local side. Selector-only faces with no mate are kept.
    static std::vector<std::shared_ptr<Sideset>> dedup_sidesets_shared_faces(
            const std::shared_ptr<Mesh>           &mesh,
            std::vector<std::shared_ptr<Sideset>>  sidesets) {
        static const int NN = LocalSideTable::MAX_NUM_NODES_PER_SIDE;
        const size_t     n_ss = sidesets.size();

        ptrdiff_t n_rec = 0;
        for (size_t i = 0; i < n_ss; ++i) {
            if (sidesets[i]) {
                n_rec += sidesets[i]->size();
            }
        }
        if (n_rec == 0) {
            return sidesets;
        }

        auto nodes   = create_host_buffer<idx_t>((size_t)n_rec * (size_t)NN);
        auto nn      = create_host_buffer<i16>((size_t)n_rec);
        auto b       = create_host_buffer<block_idx_t>((size_t)n_rec);
        auto e       = create_host_buffer<element_idx_t>((size_t)n_rec);
        auto s       = create_host_buffer<i16>((size_t)n_rec);
        auto ss_i    = create_host_buffer<i32>((size_t)n_rec);
        auto local_i = create_host_buffer<i64>((size_t)n_rec);
        auto order   = create_host_buffer<i64>((size_t)n_rec);

        idx_t         *const d_nodes   = nodes->data();
        i16           *const d_nn      = nn->data();
        block_idx_t   *const d_b       = b->data();
        element_idx_t *const d_e       = e->data();
        i16           *const d_s       = s->data();
        i32           *const d_ss_i    = ss_i->data();
        i64           *const d_local_i = local_i->data();

        ptrdiff_t w = 0;
        for (size_t i = 0; i < n_ss; ++i) {
            const auto &ss = sidesets[i];
            if (!ss || ss->size() == 0) {
                continue;
            }
            const element_idx_t *parent = ss->parent()->data();
            const i16           *lfi    = ss->lfi()->data();
            auto                 block  = mesh->block(ss->block_id());
            const enum ElemType  et     = block->element_type();
            const bool           is_ss  = is_semistructured_type(et);
            const enum ElemType  family = is_ss ? ss_source_family(et) : et;
            int                  corners[8];
            int                  n_corners = 0;
            if (is_ss) {
                ss_source_family_corners(family, semistructured_level(et), corners, &n_corners);
            }
            (void)n_corners;
            LocalSideTable lst;
            if (lst.fill(family) != SMESH_SUCCESS) {
                LocalSideTable::report_unsupported("dedup_sidesets_shared_faces", family);
                return {};
            }
            idx_t *const *const elems = block->elements()->data();
            const ptrdiff_t     n_s   = ss->size();
            const block_idx_t   bid   = ss->block_id();
            for (ptrdiff_t j = 0; j < n_s; ++j) {
                int nni = 0;
                fill_side_node_key(is_ss, corners, lst, elems, parent[j], lfi[j], d_nodes + w * NN, &nni);
                d_nn[w]      = (i16)nni;
                d_b[w]       = bid;
                d_e[w]       = parent[j];
                d_s[w]       = lfi[j];
                d_ss_i[w]    = (i32)i;
                d_local_i[w] = (i64)j;
                ++w;
            }
        }

        auto same_face = [&](const i64 a, const i64 b) -> int {
            if (d_nn[a] != d_nn[b]) {
                return 0;
            }
            const idx_t *na = d_nodes + a * NN;
            const idx_t *nb = d_nodes + b * NN;
            for (int k = 0; k < d_nn[a]; ++k) {
                if (na[k] != nb[k]) {
                    return 0;
                }
            }
            return 1;
        };

        argsort(n_rec, order->data(), [&](const i64 a, const i64 b) {
            if (d_nn[a] != d_nn[b]) {
                return d_nn[a] < d_nn[b];
            }
            const idx_t *na = d_nodes + a * NN;
            const idx_t *nb = d_nodes + b * NN;
            for (int k = 0; k < d_nn[a]; ++k) {
                if (na[k] != nb[k]) {
                    return na[k] < nb[k];
                }
            }
            if (d_b[a] != d_b[b]) {
                return d_b[a] < d_b[b];
            }
            if (d_e[a] != d_e[b]) {
                return d_e[a] < d_e[b];
            }
            return d_s[a] < d_s[b];
        });

        const i64 *const d_ord = order->data();
        auto             keep  = create_host_buffer<u8>((size_t)n_rec);
        u8 *const        d_keep = keep->data();
        std::vector<ptrdiff_t> ss_off(n_ss + 1, 0);
        for (size_t i = 0; i < n_ss; ++i) {
            ss_off[i + 1] = ss_off[i] + (sidesets[i] ? sidesets[i]->size() : 0);
        }

        for (ptrdiff_t i = 0; i < n_rec; ++i) {
            const i64 oi = d_ord[i];
            if (i == 0 || !same_face(oi, d_ord[i - 1])) {
                d_keep[ss_off[(size_t)d_ss_i[oi]] + (ptrdiff_t)d_local_i[oi]] = 1;
            }
        }

        for (size_t i = 0; i < n_ss; ++i) {
            auto &ss = sidesets[i];
            if (!ss || ss->size() == 0) {
                continue;
            }
            const ptrdiff_t      n_sides = ss->size();
            const u8            *keep_i  = d_keep + ss_off[i];
            ptrdiff_t            n_keep  = 0;
            for (ptrdiff_t j = 0; j < n_sides; ++j) {
                n_keep += keep_i[j];
            }
            if (n_keep == n_sides) {
                continue;
            }
            auto                 np   = create_host_buffer<element_idx_t>((size_t)n_keep);
            auto                 nl   = create_host_buffer<i16>((size_t)n_keep);
            const element_idx_t *p    = ss->parent()->data();
            const i16           *ls   = ss->lfi()->data();
            element_idx_t       *d_np = n_keep > 0 ? np->data() : nullptr;
            i16                 *d_nl = n_keep > 0 ? nl->data() : nullptr;
            ptrdiff_t            o    = 0;
            for (ptrdiff_t j = 0; j < n_sides; ++j) {
                if (keep_i[j]) {
                    d_np[o] = p[j];
                    d_nl[o] = ls[j];
                    ++o;
                }
            }
            ss = sideset_from_buffers(mesh, ss->block_id(), np, nl, ss->element_mapping());
        }

        return sidesets;
    }

    static int split_mixed_arity_into(const std::shared_ptr<Mesh>    &mesh,
                                      const std::shared_ptr<Sideset> &sideset,
                                      std::shared_ptr<Sideset>       *parts) {
        if (!mesh || !sideset) {
            return 0;
        }
        const enum ElemType et = mesh->element_type(sideset->block_id());
        if (elem_sides_homogeneous(et) || sideset->size() == 0) {
            parts[0] = sideset;
            return 1;
        }

        const ptrdiff_t      n      = sideset->size();
        const element_idx_t *parent = sideset->parent()->data();
        const i16           *lfi    = sideset->lfi()->data();
        ptrdiff_t            n_tri  = 0;
        ptrdiff_t            n_quad = 0;
        for (ptrdiff_t i = 0; i < n; ++i) {
            if (side_type(et, lfi[i]) == TRI3) {
                ++n_tri;
            } else {
                ++n_quad;
            }
        }

        int n_parts = 0;
        for (int want_tri = 1; want_tri >= 0; --want_tri) {
            const ptrdiff_t n_keep = want_tri ? n_tri : n_quad;
            if (n_keep == 0) {
                continue;
            }
            auto           np   = create_host_buffer<element_idx_t>((size_t)n_keep);
            auto           nl   = create_host_buffer<i16>((size_t)n_keep);
            element_idx_t *d_np = np->data();
            i16           *d_nl = nl->data();
            ptrdiff_t      o    = 0;
            for (ptrdiff_t i = 0; i < n; ++i) {
                const int is_tri = (side_type(et, lfi[i]) == TRI3);
                if (is_tri == want_tri) {
                    d_np[o] = parent[i];
                    d_nl[o] = lfi[i];
                    ++o;
                }
            }
            parts[n_parts++] =
                    sideset_from_buffers(mesh, sideset->block_id(), np, nl, sideset->element_mapping());
        }
        if (n_parts == 0) {
            parts[0] = sideset;
            return 1;
        }
        return n_parts;
    }

    std::vector<std::shared_ptr<Sideset>> split_mixed_arity_sideset(const std::shared_ptr<Mesh>    &mesh,
                                                                    const std::shared_ptr<Sideset> &sideset) {
        std::shared_ptr<Sideset> parts[2];
        const int                n = split_mixed_arity_into(mesh, sideset, parts);
        std::vector<std::shared_ptr<Sideset>> out;
        out.reserve((size_t)n);
        for (int i = 0; i < n; ++i) {
            out.push_back(parts[i]);
        }
        return out;
    }

    static std::shared_ptr<Buffer<idx_t *>> concat_surface_soa(
            const std::vector<std::shared_ptr<Buffer<idx_t *>>> &parts) {
        const int n_parts = (int)parts.size();
        if (n_parts <= 0) {
            return nullptr;
        }
        const int nnxs = (int)parts[0]->extent(0);
        ptrdiff_t n    = 0;
        for (int i = 0; i < n_parts; ++i) {
            if (!parts[i] || (int)parts[i]->extent(0) != nnxs) {
                return nullptr;
            }
            n += (ptrdiff_t)parts[i]->extent(1);
        }
        auto      out   = create_host_buffer<idx_t>(nnxs, n);
        idx_t   **d_out = out->data();
        ptrdiff_t off   = 0;
        for (int i = 0; i < n_parts; ++i) {
            const ptrdiff_t  ne     = (ptrdiff_t)parts[i]->extent(1);
            const idx_t *const *d_in = parts[i]->data();
            for (int d = 0; d < nnxs; ++d) {
                if (ne > 0) {
                    memcpy(d_out[d] + off, d_in[d], (size_t)ne * sizeof(idx_t));
                }
            }
            off += ne;
        }
        return out;
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
            if (!block_name_selected(block_names, block->name())) {
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
            if (lst.fill(family) != SMESH_SUCCESS) {
                LocalSideTable::report_unsupported("Sideset::create_from_selector", family);
                return {};
            }

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

        return dedup_sidesets_shared_faces(mesh, filter_sidesets_to_owned(mesh, sidesets));
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
            if (!block_name_selected(block_names, block->name())) {
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
            if (lst.fill(family) != SMESH_SUCCESS) {
                LocalSideTable::report_unsupported("Sideset::create_from_batch_selector", family);
                return {};
            }

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

        return dedup_sidesets_shared_faces(mesh, filter_sidesets_to_owned(mesh, sidesets));
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

        bool        has_size     = false;
        bool        has_block_id = false;
        ptrdiff_t   meta_size    = 0;
        block_idx_t meta_block_id = 0;
        if (parse_sideset_meta(folder, &has_size, &meta_size, &has_block_id, &meta_block_id) != SMESH_SUCCESS) {
            return SMESH_FAILURE;
        }
        if (has_size && meta_size != nlocal) {
            SMESH_ERROR("Sideset meta size %td does not match array length %td at %s\n",
                        meta_size,
                        nlocal,
                        folder.c_str());
            ret = SMESH_FAILURE;
        }

        impl_->block_id = has_block_id ? meta_block_id : block_id;
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
                    n_sides  = static_cast<i64>(this->size());
                    block_id = impl_->block_id;
                }
            }
            mesh->comm()->broadcast(&n_sides, 1, 0);
            if (n_sides < 0) {
                return SMESH_FAILURE;
            }

            mesh->comm()->broadcast(&block_id, 1, 0);
            impl_->block_id = block_id;
            impl_->comm     = mesh->comm();

            auto parent_buf = create_host_buffer<element_idx_t>((size_t)n_sides);
            auto lfi_buf    = create_host_buffer<i16>((size_t)n_sides);
            if (rank == 0 && n_sides > 0) {
                std::memcpy(parent_buf->data(),
                            impl_->parent->data(),
                            (size_t)n_sides * sizeof(element_idx_t));
                std::memcpy(lfi_buf->data(), impl_->lfi->data(), (size_t)n_sides * sizeof(i16));
            }
            if (n_sides > 0) {
                mesh->comm()->broadcast(parent_buf->data(), (int)n_sides, 0);
                mesh->comm()->broadcast(lfi_buf->data(), (int)n_sides, 0);
            }

            impl_->parent = parent_buf;
            impl_->lfi    = lfi_buf;
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
        const large_idx_t *d_emap = (n_owned > 0) ? elem_map->data() : nullptr;
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
            impl_->lfi             = create_host_buffer<i16>(0);
            impl_->element_mapping = block->element_mapping();
            return SMESH_SUCCESS;
        }

        auto            invert = create_host_buffer<element_idx_t>((size_t)n_global);
        element_idx_t  *d_inv  = invert->data();
        for (large_idx_t i = 0; i < n_global; ++i) {
            d_inv[i] = invalid_idx<element_idx_t>();
        }
        for (ptrdiff_t i = 0; i < n_owned; ++i) {
            const large_idx_t serial_id = d_emap[i];
            if (serial_id >= 0 && serial_id < n_global) {
                d_inv[serial_id] = (element_idx_t)i;
            }
        }

        const ptrdiff_t      n_sides = size();
        const element_idx_t *parent  = impl_->parent->data();
        const i16           *lfi     = impl_->lfi->data();

        ptrdiff_t n_keep = 0;
        for (ptrdiff_t i = 0; i < n_sides; ++i) {
            const large_idx_t serial_parent = (large_idx_t)parent[i];
            if (serial_parent >= 0 && serial_parent < n_global && d_inv[serial_parent] != invalid_idx<element_idx_t>()) {
                ++n_keep;
            }
        }

        auto new_parent = create_host_buffer<element_idx_t>((size_t)n_keep);
        auto new_lfi    = create_host_buffer<i16>((size_t)n_keep);
        element_idx_t *d_np = n_keep > 0 ? new_parent->data() : nullptr;
        i16           *d_nl = n_keep > 0 ? new_lfi->data() : nullptr;
        ptrdiff_t w     = 0;
        for (ptrdiff_t i = 0; i < n_sides; ++i) {
            const large_idx_t serial_parent = (large_idx_t)parent[i];
            if (serial_parent >= 0 && serial_parent < n_global) {
                const element_idx_t local_e = d_inv[serial_parent];
                if (local_e != invalid_idx<element_idx_t>()) {
                    d_np[w] = local_e;
                    d_nl[w] = lfi[i];
                    ++w;
                }
            }
        }

        impl_->parent          = new_parent;
        impl_->lfi             = new_lfi;
        impl_->element_mapping = block->element_mapping();

        return SMESH_SUCCESS;
#endif

        SMESH_ERROR("Sideset::redistribute is not supported for distributed runs\n");
        return SMESH_FAILURE;
    }

    std::shared_ptr<Buffer<element_idx_t>> Sideset::parent() const { return impl_->parent; }
    std::shared_ptr<Buffer<i16>>           Sideset::lfi() const { return impl_->lfi; }
    block_idx_t                            Sideset::block_id() const { return impl_->block_id; }
    SharedBuffer<large_idx_t>              Sideset::element_mapping() const { return impl_->element_mapping; }

    int Sideset::remap_parents(const element_idx_t *old_to_new, ptrdiff_t n) {
        if (!old_to_new || n < 0 || !impl_->parent) {
            SMESH_ERROR("Sideset::remap_parents: invalid arguments\n");
            return SMESH_FAILURE;
        }
        if (impl_->parent->mem_space() != MEMORY_SPACE_HOST) {
            SMESH_ERROR("Sideset::remap_parents requires host parent buffer\n");
            return SMESH_FAILURE;
        }
        auto            p   = impl_->parent->data();
        const ptrdiff_t nss = size();
        for (ptrdiff_t i = 0; i < nss; ++i) {
            const element_idx_t old = p[i];
            if (old < 0 || static_cast<ptrdiff_t>(old) >= n) {
                SMESH_ERROR("Sideset::remap_parents: parent %td out of range [0, %td)\n",
                            static_cast<ptrdiff_t>(old),
                            n);
                return SMESH_FAILURE;
            }
            p[i] = old_to_new[old];
        }
        return SMESH_SUCCESS;
    }

    int Sideset::apply_element_permutation(const element_idx_t *new_to_old, ptrdiff_t n) {
        if (!new_to_old || n < 0) {
            SMESH_ERROR("Sideset::apply_element_permutation: invalid arguments\n");
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
                SMESH_ERROR("Sideset::apply_element_permutation: gather index %td out of range\n",
                            (ptrdiff_t)old);
                return SMESH_FAILURE;
            }
            d_otn[old] = (element_idx_t)neu;
        }
        for (ptrdiff_t old = 0; old < n; ++old) {
            if (d_otn[old] == invalid_idx<element_idx_t>()) {
                SMESH_ERROR("Sideset::apply_element_permutation: permutation is not a bijection\n");
                return SMESH_FAILURE;
            }
        }
        return remap_parents(d_otn, n);
    }

    int remap_sidesets(const std::vector<std::shared_ptr<Sideset>> &sidesets,
                       block_idx_t                                  block_id,
                       const element_idx_t                         *old_to_new,
                       ptrdiff_t                                    n) {
        for (auto &ss : sidesets) {
            if (!ss || ss->block_id() != block_id) {
                continue;
            }
            if (ss->remap_parents(old_to_new, n) != SMESH_SUCCESS) {
                return SMESH_FAILURE;
            }
        }
        return SMESH_SUCCESS;
    }

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
            std::vector<std::shared_ptr<Buffer<idx_t>>> parts;
            parts.reserve(sidesets.size());
            ptrdiff_t n_ids = 0;
            for (size_t k = 0; k < sidesets.size(); ++k) {
                auto ns = create_nodeset_from_sideset(mesh, sidesets[k]);
                if (ns && ns->size() > 0) {
                    n_ids += (ptrdiff_t)ns->size();
                    parts.push_back(std::move(ns));
                }
            }
            auto ids = create_host_buffer<idx_t>((size_t)n_ids);
            idx_t    *d_ids = n_ids > 0 ? ids->data() : nullptr;
            ptrdiff_t w     = 0;
            for (size_t k = 0; k < parts.size(); ++k) {
                memcpy(d_ids + w, parts[k]->data(), parts[k]->size() * sizeof(idx_t));
                w += (ptrdiff_t)parts[k]->size();
            }
            const ptrdiff_t n = w > 0 ? (ptrdiff_t)sort_and_unique(d_ids, (size_t)w) : 0;
            return n > 0 ? view(ids, 0, (size_t)n) : ids;
        }

        const block_idx_t                 n_sidesets = (block_idx_t)sidesets.size();
        std::vector<enum ElemType>        element_type((size_t)n_sidesets);
        std::vector<idx_t **>             elems((size_t)n_sidesets);
        std::vector<ptrdiff_t>            n_surf_elements((size_t)n_sidesets);
        std::vector<element_idx_t *>      parent_element((size_t)n_sidesets);
        std::vector<i16 *>                side_idx((size_t)n_sidesets);

        for (block_idx_t k = 0; k < n_sidesets; k++) {
            auto ss    = sidesets[k];
            auto block = mesh->block(ss->block_id());

            element_type[k]    = block->element_type();
            elems[k]           = block->elements()->data();
            n_surf_elements[k] = ss->parent()->size();
            parent_element[k]  = ss->parent()->data();
            side_idx[k]        = ss->lfi()->data();
        }

        ptrdiff_t n_nodes_out = 0;
        idx_t    *nodes_out   = nullptr;

        extract_nodeset_from_sidesets(n_sidesets,
                                      element_type.data(),
                                      elems.data(),
                                      n_surf_elements.data(),
                                      parent_element.data(),
                                      side_idx.data(),
                                      &n_nodes_out,
                                      &nodes_out);

        return manage_host_buffer(n_nodes_out, nodes_out);
    }

    static int hex_face_child_on(const int L, const int lfi, const int xi, const int yi, const int zi) {
        switch (lfi) {
            case 0:
                return yi == 0;
            case 1:
                return xi == L - 1;
            case 2:
                return yi == L - 1;
            case 3:
                return xi == 0;
            case 4:
                return zi == 0;
            case 5:
                return zi == L - 1;
            default:
                return 0;
        }
    }

    static int quad_edge_child_on(const int L, const int lfi, const int xi, const int yi) {
        switch (lfi) {
            case 0:
                return yi == 0;
            case 1:
                return xi == L - 1;
            case 2:
                return yi == L - 1;
            case 3:
                return xi == 0;
            default:
                return 0;
        }
    }

    /// WEDGE6 faces: 0 y=0 QUAD, 1 hypotenuse QUAD, 2 x=0 QUAD, 3 z=0 TRI, 4 z=L TRI.
    /// `down` is the opposite-orientation microtriangle in the (x,y) plane.
    static int wedge_child_on(const int L, const int lfi, const int xi, const int yi, const int zi,
                              const int down) {
        switch (lfi) {
            case 0:
                return !down && yi == 0;
            case 1:
                return !down && xi + yi == L - 1;
            case 2:
                return !down && xi == 0;
            case 3:
                return zi == 0;
            case 4:
                return zi == L - 1;
            default:
                return 0;
        }
    }

    std::shared_ptr<Sideset> map_sideset_through_refine(const std::shared_ptr<Mesh>    &coarse_mesh,
                                                        const std::shared_ptr<Sideset> &coarse_ss,
                                                        const std::shared_ptr<Mesh>    &fine_mesh) {
        if (!coarse_mesh || !coarse_ss || !fine_mesh) {
            SMESH_ERROR("map_sideset_through_refine: null argument\n");
            return nullptr;
        }
        const block_idx_t bid          = coarse_ss->block_id();
        auto              coarse_block = coarse_mesh->block(bid);
        auto              fine_block   = fine_mesh->block(bid);
        if (!coarse_block || !fine_block) {
            SMESH_ERROR("map_sideset_through_refine: missing block %d\n", static_cast<int>(bid));
            return nullptr;
        }

        const enum ElemType coarse_type = coarse_block->element_type();
        const enum ElemType fine_type   = fine_block->element_type();
        const ptrdiff_t     n_coarse    = coarse_block->n_elements();
        const ptrdiff_t     n_fine      = fine_block->n_elements();
        if (n_coarse <= 0) {
            SMESH_ERROR("map_sideset_through_refine: empty coarse block\n");
            return nullptr;
        }

        ptrdiff_t factor = 0;
        if (n_fine % n_coarse == 0) {
            factor = n_fine / n_coarse;
        }

        auto unsupported = [&]() -> std::shared_ptr<Sideset> {
            fprintf(stderr,
                    "map_sideset_through_refine: unsupported for %s -> %s "
                    "(n_coarse=%td n_fine=%td factor=%td)\n",
                    type_to_string(coarse_type),
                    type_to_string(fine_type),
                    n_coarse,
                    n_fine,
                    factor);
            return nullptr;
        };

        if (is_semistructured_type(coarse_type) || is_semistructured_type(fine_type)) {
            fprintf(stderr,
                    "map_sideset_through_refine: SS meshes keep (parent, lfi); do not remap (%s -> %s)\n",
                    type_to_string(coarse_type),
                    type_to_string(fine_type));
            return nullptr;
        }

        if (factor == 1) {
            if (elem_num_sides(coarse_type) != elem_num_sides(fine_type)) {
                return unsupported();
            }
            auto parent = create_host_buffer<element_idx_t>(coarse_ss->size());
            auto lfi    = create_host_buffer<i16>(coarse_ss->size());
            std::memcpy(parent->data(),
                        coarse_ss->parent()->data(),
                        static_cast<size_t>(coarse_ss->size()) * sizeof(element_idx_t));
            std::memcpy(lfi->data(), coarse_ss->lfi()->data(), static_cast<size_t>(coarse_ss->size()) * sizeof(i16));
            return Sideset::create(fine_mesh->comm(),
                                   parent,
                                   lfi,
                                   bid,
                                   fine_mesh->comm()->size() > 1 ? fine_block->element_mapping() : nullptr);
        }

        const ptrdiff_t      n_ss = coarse_ss->size();
        const element_idx_t *cp   = coarse_ss->parent()->data();
        const i16           *cl   = coarse_ss->lfi()->data();
        SharedBuffer<element_idx_t> out_p;
        SharedBuffer<i16>           out_l;
        ptrdiff_t                   n_out = 0;

        if (coarse_type == HEX8 && fine_type == HEX8) {
            if (factor <= 0) {
                return unsupported();
            }
            int L = 1;
            while ((ptrdiff_t)L * (ptrdiff_t)L * (ptrdiff_t)L < factor) {
                ++L;
            }
            if ((ptrdiff_t)L * (ptrdiff_t)L * (ptrdiff_t)L != factor) {
                return unsupported();
            }
            n_out = n_ss * (ptrdiff_t)L * (ptrdiff_t)L;
            out_p = create_host_buffer<element_idx_t>((size_t)n_out);
            out_l = create_host_buffer<i16>((size_t)n_out);
            element_idx_t *d_out_p = n_out > 0 ? out_p->data() : nullptr;
            i16           *d_out_l = n_out > 0 ? out_l->data() : nullptr;
            ptrdiff_t w = 0;
            for (ptrdiff_t i = 0; i < n_ss; ++i) {
                const element_idx_t e   = cp[i];
                const i16           lfi = cl[i];
                if (e < 0 || (ptrdiff_t)e >= n_coarse) {
                    SMESH_ERROR("map_sideset_through_refine: parent out of range\n");
                    return nullptr;
                }
                if (lfi < 0 || lfi > 5) {
                    SMESH_ERROR("map_sideset_through_refine: HEX8 side %d is invalid\n", (int)lfi);
                    return nullptr;
                }
                for (int zi = 0; zi < L; ++zi) {
                    for (int yi = 0; yi < L; ++yi) {
                        for (int xi = 0; xi < L; ++xi) {
                            if (hex_face_child_on(L, lfi, xi, yi, zi)) {
                                d_out_p[w] = (element_idx_t)((ptrdiff_t)e * factor + zi * L * L + yi * L + xi);
                                d_out_l[w] = lfi;
                                ++w;
                            }
                        }
                    }
                }
            }
            n_out = w;
        } else if ((coarse_type == QUAD4 || coarse_type == QUADSHELL4) &&
                   (fine_type == QUAD4 || fine_type == QUADSHELL4)) {
            if (factor <= 0) {
                return unsupported();
            }
            int L = 1;
            while ((ptrdiff_t)L * (ptrdiff_t)L < factor) {
                ++L;
            }
            if ((ptrdiff_t)L * (ptrdiff_t)L != factor) {
                return unsupported();
            }
            n_out = n_ss * (ptrdiff_t)L;
            out_p = create_host_buffer<element_idx_t>((size_t)n_out);
            out_l = create_host_buffer<i16>((size_t)n_out);
            element_idx_t *d_out_p = n_out > 0 ? out_p->data() : nullptr;
            i16           *d_out_l = n_out > 0 ? out_l->data() : nullptr;
            ptrdiff_t w = 0;
            for (ptrdiff_t i = 0; i < n_ss; ++i) {
                const element_idx_t e   = cp[i];
                const i16           lfi = cl[i];
                if (e < 0 || (ptrdiff_t)e >= n_coarse) {
                    SMESH_ERROR("map_sideset_through_refine: parent out of range\n");
                    return nullptr;
                }
                if (lfi < 0 || lfi > 3) {
                    SMESH_ERROR("map_sideset_through_refine: QUAD side %d is invalid\n", (int)lfi);
                    return nullptr;
                }
                for (int yi = 0; yi < L; ++yi) {
                    for (int xi = 0; xi < L; ++xi) {
                        if (quad_edge_child_on(L, lfi, xi, yi)) {
                            d_out_p[w] = (element_idx_t)((ptrdiff_t)e * factor + yi * L + xi);
                            d_out_l[w] = lfi;
                            ++w;
                        }
                    }
                }
            }
            n_out = w;
        } else if (coarse_type == WEDGE6 && fine_type == WEDGE6) {
            if (factor <= 0) {
                return unsupported();
            }
            int L = 1;
            while ((ptrdiff_t)L * (ptrdiff_t)L * (ptrdiff_t)L < factor) {
                ++L;
            }
            if ((ptrdiff_t)L * (ptrdiff_t)L * (ptrdiff_t)L != factor) {
                return unsupported();
            }
            n_out = n_ss * (ptrdiff_t)L * (ptrdiff_t)L;
            out_p = create_host_buffer<element_idx_t>((size_t)n_out);
            out_l = create_host_buffer<i16>((size_t)n_out);
            element_idx_t *d_out_p = n_out > 0 ? out_p->data() : nullptr;
            i16           *d_out_l = n_out > 0 ? out_l->data() : nullptr;
            ptrdiff_t w = 0;
            for (ptrdiff_t i = 0; i < n_ss; ++i) {
                const element_idx_t e   = cp[i];
                const i16           lfi = cl[i];
                if (e < 0 || (ptrdiff_t)e >= n_coarse) {
                    SMESH_ERROR("map_sideset_through_refine: parent out of range\n");
                    return nullptr;
                }
                if (lfi < 0 || lfi > 4) {
                    SMESH_ERROR("map_sideset_through_refine: WEDGE6 side %d is invalid\n", (int)lfi);
                    return nullptr;
                }
                int le = 0;
                for (int zi = 0; zi < L; ++zi) {
                    for (int yi = 0; yi < L; ++yi) {
                        for (int xi = 0; xi < L - yi; ++xi) {
                            if (wedge_child_on(L, lfi, xi, yi, zi, 0)) {
                                d_out_p[w] = (element_idx_t)((ptrdiff_t)e * factor + le);
                                d_out_l[w] = lfi;
                                ++w;
                            }
                            ++le;
                            if (xi + yi + 1 < L) {
                                if (wedge_child_on(L, lfi, xi, yi, zi, 1)) {
                                    d_out_p[w] = (element_idx_t)((ptrdiff_t)e * factor + le);
                                    d_out_l[w] = lfi;
                                    ++w;
                                }
                                ++le;
                            }
                        }
                    }
                }
            }
            n_out = w;
        } else if (coarse_type == TET4 && fine_type == TET4) {
            if (factor <= 0) {
                return unsupported();
            }
            ptrdiff_t f  = 1;
            int       lv = 0;
            while (f < factor) {
                f *= 8;
                ++lv;
            }
            if (f != factor) {
                return unsupported();
            }
            static const int tet_face_child[4][4] = {{0, 1, 3, 5}, {1, 2, 3, 6}, {0, 2, 3, 7}, {0, 1, 2, 4}};
            ptrdiff_t        expand               = 1;
            for (int s = 0; s < lv; ++s) {
                expand *= 4;
            }
            n_out = n_ss * expand;
            out_p = create_host_buffer<element_idx_t>((size_t)n_out);
            out_l = create_host_buffer<i16>((size_t)n_out);
            auto cur_p = create_host_buffer<element_idx_t>((size_t)expand);
            auto cur_l = create_host_buffer<i16>((size_t)expand);
            auto nxt_p = create_host_buffer<element_idx_t>((size_t)expand);
            auto nxt_l = create_host_buffer<i16>((size_t)expand);
            element_idx_t *d_out_p = n_out > 0 ? out_p->data() : nullptr;
            i16           *d_out_l = n_out > 0 ? out_l->data() : nullptr;
            element_idx_t *d_cur_p = cur_p->data();
            i16           *d_cur_l = cur_l->data();
            element_idx_t *d_nxt_p = nxt_p->data();
            i16           *d_nxt_l = nxt_l->data();
            ptrdiff_t      w       = 0;
            for (ptrdiff_t i = 0; i < n_ss; ++i) {
                d_cur_p[0]      = cp[i];
                d_cur_l[0]      = cl[i];
                ptrdiff_t n_cur = 1;
                for (int step = 0; step < lv; ++step) {
                    ptrdiff_t n_nxt = 0;
                    for (ptrdiff_t k = 0; k < n_cur; ++k) {
                        const i16 s = d_cur_l[k];
                        if (s < 0 || s > 3) {
                            SMESH_ERROR("map_sideset_through_refine: TET4 side %d is invalid\n", (int)s);
                            return nullptr;
                        }
                        for (int c = 0; c < 4; ++c) {
                            d_nxt_p[n_nxt] = (element_idx_t)((ptrdiff_t)d_cur_p[k] * 8 + tet_face_child[s][c]);
                            d_nxt_l[n_nxt] = s;
                            ++n_nxt;
                        }
                    }
                    element_idx_t *tp = d_cur_p;
                    i16           *tl = d_cur_l;
                    d_cur_p           = d_nxt_p;
                    d_cur_l           = d_nxt_l;
                    d_nxt_p           = tp;
                    d_nxt_l           = tl;
                    n_cur             = n_nxt;
                }
                memcpy(d_out_p + w, d_cur_p, (size_t)n_cur * sizeof(element_idx_t));
                memcpy(d_out_l + w, d_cur_l, (size_t)n_cur * sizeof(i16));
                w += n_cur;
            }
            n_out = w;
        } else if ((coarse_type == TRI3 || coarse_type == TRISHELL3) && fine_type == coarse_type) {
            if (factor <= 0) {
                return unsupported();
            }
            ptrdiff_t f  = 1;
            int       lv = 0;
            while (f < factor) {
                f *= 4;
                ++lv;
            }
            if (f != factor) {
                return unsupported();
            }
            static const int tri_edge_child[3][2] = {{0, 1}, {1, 2}, {2, 0}};
            ptrdiff_t        expand               = 1;
            for (int s = 0; s < lv; ++s) {
                expand *= 2;
            }
            n_out = n_ss * expand;
            out_p = create_host_buffer<element_idx_t>((size_t)n_out);
            out_l = create_host_buffer<i16>((size_t)n_out);
            auto cur_p = create_host_buffer<element_idx_t>((size_t)expand);
            auto cur_l = create_host_buffer<i16>((size_t)expand);
            auto nxt_p = create_host_buffer<element_idx_t>((size_t)expand);
            auto nxt_l = create_host_buffer<i16>((size_t)expand);
            element_idx_t *d_out_p = n_out > 0 ? out_p->data() : nullptr;
            i16           *d_out_l = n_out > 0 ? out_l->data() : nullptr;
            element_idx_t *d_cur_p = cur_p->data();
            i16           *d_cur_l = cur_l->data();
            element_idx_t *d_nxt_p = nxt_p->data();
            i16           *d_nxt_l = nxt_l->data();
            ptrdiff_t      w       = 0;
            for (ptrdiff_t i = 0; i < n_ss; ++i) {
                d_cur_p[0]      = cp[i];
                d_cur_l[0]      = cl[i];
                ptrdiff_t n_cur = 1;
                for (int step = 0; step < lv; ++step) {
                    ptrdiff_t n_nxt = 0;
                    for (ptrdiff_t k = 0; k < n_cur; ++k) {
                        const i16 s = d_cur_l[k];
                        if (s < 0 || s > 2) {
                            SMESH_ERROR("map_sideset_through_refine: TRI side %d is invalid\n", (int)s);
                            return nullptr;
                        }
                        for (int c = 0; c < 2; ++c) {
                            d_nxt_p[n_nxt] = (element_idx_t)((ptrdiff_t)d_cur_p[k] * 4 + tri_edge_child[s][c]);
                            d_nxt_l[n_nxt] = s;
                            ++n_nxt;
                        }
                    }
                    element_idx_t *tp = d_cur_p;
                    i16           *tl = d_cur_l;
                    d_cur_p           = d_nxt_p;
                    d_cur_l           = d_nxt_l;
                    d_nxt_p           = tp;
                    d_nxt_l           = tl;
                    n_cur             = n_nxt;
                }
                memcpy(d_out_p + w, d_cur_p, (size_t)n_cur * sizeof(element_idx_t));
                memcpy(d_out_l + w, d_cur_l, (size_t)n_cur * sizeof(i16));
                w += n_cur;
            }
            n_out = w;
        } else if ((coarse_type == EDGE2 || coarse_type == EDGESHELL2) && fine_type == coarse_type) {
            if (factor <= 0) {
                return unsupported();
            }
            ptrdiff_t f = 1;
            while (f < factor) {
                f *= 2;
            }
            if (f != factor) {
                return unsupported();
            }
            n_out = n_ss;
            out_p = create_host_buffer<element_idx_t>((size_t)n_out);
            out_l = create_host_buffer<i16>((size_t)n_out);
            element_idx_t *d_out_p = n_out > 0 ? out_p->data() : nullptr;
            i16           *d_out_l = n_out > 0 ? out_l->data() : nullptr;
            for (ptrdiff_t i = 0; i < n_ss; ++i) {
                const element_idx_t e   = cp[i];
                const i16           lfi = cl[i];
                if (e < 0 || (ptrdiff_t)e >= n_coarse) {
                    SMESH_ERROR("map_sideset_through_refine: parent out of range\n");
                    return nullptr;
                }
                if (lfi != 0 && lfi != 1) {
                    SMESH_ERROR("map_sideset_through_refine: EDGE side %d is invalid\n", (int)lfi);
                    return nullptr;
                }
                const ptrdiff_t child =
                        (ptrdiff_t)e * factor + ((lfi == 0) ? (ptrdiff_t)0 : (factor - 1));
                d_out_p[i] = (element_idx_t)child;
                d_out_l[i] = lfi;
            }
        } else {
            return unsupported();
        }

        if (n_out > 0 && (ptrdiff_t)out_p->size() != n_out) {
            out_p = view(out_p, 0, (size_t)n_out);
            out_l = view(out_l, 0, (size_t)n_out);
        }
        return Sideset::create(fine_mesh->comm(),
                               out_p,
                               out_l,
                               bid,
                               fine_mesh->comm()->size() > 1 ? fine_block->element_mapping() : nullptr);
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
            const i16 *d_lfi        = sideset->lfi()->data();
            const int  first        = d_lfi[0];
            const bool is_quad_side = (family == WEDGE6) ? (first < 3) : (first >= 4);
            for (ptrdiff_t i = 0; i < n_surf; ++i) {
                const int  s = d_lfi[i];
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
        if (!LocalSideTable::supported(block->element_type())) {
            LocalSideTable::report_unsupported("create_surface_from_sideset", block->element_type());
            return {INVALID, nullptr};
        }
        const ptrdiff_t n_surf = sideset->parent()->size();
        if (n_surf > 0 && !elem_sides_homogeneous(block->element_type())) {
            const i16          *d_lfi = sideset->lfi()->data();
            const enum ElemType first = side_type(block->element_type(), d_lfi[0]);
            for (ptrdiff_t i = 1; i < n_surf; ++i) {
                if (side_type(block->element_type(), d_lfi[i]) != first) {
                    fprintf(stderr,
                            "create_surface_from_sideset: mixed TRI/QUAD sideset; "
                            "use create_surfaces_from_sidesets or split_mixed_arity_sideset\n");
                    return {INVALID, nullptr};
                }
            }
        }
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

        const large_idx_t *const nmap          = nmap_buf->data();
        const ptrdiff_t          n_local_nodes = mesh->n_nodes();
        const ptrdiff_t          n_local_ids   = (local && local->size() > 0) ? (ptrdiff_t)local->size() : 0;

        auto      send  = create_host_buffer<large_idx_t>((size_t)n_local_ids);
        ptrdiff_t nsend = 0;
        if (n_local_ids > 0) {
            const idx_t    *d     = local->data();
            large_idx_t    *d_send = send->data();
            for (ptrdiff_t i = 0; i < n_local_ids; ++i) {
                const idx_t n = d[i];
                if (n >= 0 && (ptrdiff_t)n < n_local_nodes) {
                    d_send[nsend++] = nmap[n];
                }
            }
            nsend = (ptrdiff_t)sort_and_unique(d_send, (size_t)nsend);
        }

        MPI_Comm         comm   = mesh->comm()->get();
        const int        nsend_i = (int)nsend;
        const int        nranks  = mesh->comm()->size();
        std::vector<int> counts((size_t)nranks, 0);
        std::vector<int> displs((size_t)nranks, 0);
        SMESH_MPI_CATCH(MPI_Allgather(&nsend_i, 1, MPI_INT, counts.data(), 1, MPI_INT, comm));
        int ntotal = 0;
        for (int r = 0; r < nranks; ++r) {
            displs[(size_t)r] = ntotal;
            ntotal += counts[(size_t)r];
        }

        auto recv = create_host_buffer<large_idx_t>((size_t)ntotal);
        SMESH_MPI_CATCH(MPI_Allgatherv(nsend > 0 ? send->data() : nullptr,
                                       nsend_i,
                                       mpi_type<large_idx_t>(),
                                       ntotal > 0 ? recv->data() : nullptr,
                                       counts.data(),
                                       displs.data(),
                                       mpi_type<large_idx_t>(),
                                       comm));

        if (ntotal <= 0 || n_local_nodes <= 0) {
            return create_host_buffer<idx_t>(0);
        }

        auto order = create_host_buffer<idx_t>((size_t)n_local_nodes);
        argsort(n_local_nodes, nmap, order->data());
        const idx_t *const d_ord = order->data();

        auto      out  = create_host_buffer<idx_t>((size_t)ntotal);
        ptrdiff_t nout = 0;
        const large_idx_t *d_recv = recv->data();
        idx_t             *d_out  = out->data();
        for (int i = 0; i < ntotal; ++i) {
            const large_idx_t gid = d_recv[i];
            const ptrdiff_t   p   = argsort_lower_bound(n_local_nodes, nmap, d_ord, gid);
            if (p < n_local_nodes && nmap[d_ord[p]] == gid) {
                d_out[nout++] = d_ord[p];
            }
        }
        nout = nout > 0 ? (ptrdiff_t)sort_and_unique(d_out, (size_t)nout) : 0;
        return nout > 0 ? view(out, 0, (size_t)nout) : create_host_buffer<idx_t>(0);
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
            if (!LocalSideTable::supported(block->element_type())) {
                LocalSideTable::report_unsupported("create_nodeset_from_sideset", block->element_type());
                return nullptr;
            }

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

    std::vector<std::pair<enum ElemType, std::shared_ptr<Buffer<idx_t *>>>> create_surfaces_from_sidesets(
            const std::shared_ptr<Mesh>                 &mesh,
            const std::vector<std::shared_ptr<Sideset>> &sidesets) {
        std::vector<enum ElemType>                                 types;
        std::vector<std::vector<std::shared_ptr<Buffer<idx_t *>>>> groups;
        for (size_t si = 0; si < sidesets.size(); ++si) {
            const auto &ss = sidesets[si];
            if (!ss || ss->size() == 0) {
                continue;
            }
            auto parts = split_mixed_arity_sideset(mesh, ss);
            for (size_t pi = 0; pi < parts.size(); ++pi) {
                const auto &part = parts[pi];
                if (!part || part->size() == 0) {
                    continue;
                }
                auto [st, surface] = create_surface_from_sideset(mesh, part);
                if (st == INVALID || !surface) {
                    return {};
                }
                size_t g = types.size();
                for (size_t t = 0; t < types.size(); ++t) {
                    if (types[t] == st) {
                        g = t;
                        break;
                    }
                }
                if (g == types.size()) {
                    types.push_back(st);
                    groups.emplace_back();
                }
                groups[g].push_back(std::move(surface));
            }
        }

        std::vector<std::pair<enum ElemType, std::shared_ptr<Buffer<idx_t *>>>> out;
        out.reserve(types.size());
        for (size_t t = 0; t < types.size(); ++t) {
            auto merged = concat_surface_soa(groups[t]);
            if (!merged) {
                return {};
            }
            out.emplace_back(types[t], std::move(merged));
        }
        return out;
    }

    std::pair<enum ElemType, std::shared_ptr<Buffer<idx_t *>>> create_surface_from_sidesets(
            const std::shared_ptr<Mesh>                 &mesh,
            const std::vector<std::shared_ptr<Sideset>> &sidesets) {
        if (sidesets.empty()) {
            return {INVALID, nullptr};
        }

        auto surfaces = create_surfaces_from_sidesets(mesh, sidesets);
        if (surfaces.size() == 1) {
            return surfaces[0];
        }
        if (surfaces.empty()) {
            return {INVALID, nullptr};
        }
        fprintf(stderr,
                "create_surface_from_sidesets: %zu distinct shell types; use create_surfaces_from_sidesets\n",
                surfaces.size());
        return {INVALID, nullptr};
    }

    void Sideset::print(std::ostream &os) const {
        os << "Sideset: " << block_id() << "\n";
        os << "Size: " << parent()->size() << "\n";

        parent()->print(os);
        lfi()->print(os);
    }

}  // namespace smesh


