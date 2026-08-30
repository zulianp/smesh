#include "smesh_restrict.hpp"
#include "smesh_buffer.hpp"
#include "smesh_device_buffer.hpp"
#include "smesh_exchange.hpp"
#include "smesh_mesh.hpp"
#include "smesh_restriction.hpp"
#include "smesh_semistructured.hpp"
#include "smesh_sshex8.hpp"
#include "smesh_sshex8_restriction.hpp"
#include "smesh_ssquad4.hpp"
#include "smesh_ssquad4_restriction.hpp"
#include "smesh_sstet4_restriction.hpp"
#include "smesh_tracer.hpp"

#include <algorithm>
#include <cstring>
#include <vector>

#ifdef SMESH_ENABLE_CUDA
#include "smesh_sshex8_restriction.cuh"
#include "smesh_ssquad4_restriction.cuh"
#include "smesh_tet4_prolongation_restriction.cuh"
#endif

namespace smesh {

    static ElemType transfer_source_family(const ElemType type) {
        return is_semistructured_type(type) ? ss_source_family(type) : type;
    }

    static ElemType ss_block_family(const Mesh &mesh, const size_t b) {
        return ss_source_family(mesh.element_type(static_cast<block_idx_t>(b)));
    }

    static bool ss_all_hex_family(const Mesh &mesh) {
        for (size_t b = 0; b < mesh.n_blocks(); ++b) {
            if (ss_block_family(mesh, b) != HEX8) {
                return false;
            }
        }
        return true;
    }

    static bool ss_has_hex_and_tet(const Mesh &mesh) {
        bool has_hex = false;
        bool has_tet = false;
        for (size_t b = 0; b < mesh.n_blocks(); ++b) {
            const auto fam = ss_block_family(mesh, b);
            has_hex |= fam == HEX8;
            has_tet |= fam == TET4;
        }
        return has_hex && has_tet;
    }

    static void add_block_incidence(const Mesh::Block &block, const ptrdiff_t nelements, uint16_t *const buff) {
        if (nelements == 0) {
            return;
        }
        const int     nxe      = block.n_nodes_per_element();
        idx_t **const elements = block.elements()->data();
        for (int d = 0; d < nxe; d++) {
            for (ptrdiff_t i = 0; i < nelements; ++i) {
                buff[elements[d][i]]++;
            }
        }
    }

    static void mark_family_nodes(const Mesh &mesh, const ElemType family, std::vector<char> &mask) {
        mask.assign((size_t)mesh.n_nodes(), 0);
        for (size_t b = 0; b < mesh.n_blocks(); ++b) {
            if (ss_block_family(mesh, b) != family) {
                continue;
            }
            auto            block = mesh.block(b);
            const ptrdiff_t ne    = block->n_elements();
            if (ne == 0) {
                continue;
            }
            const int     nxe      = block->n_nodes_per_element();
            idx_t **const elements = block->elements()->data();
            for (int d = 0; d < nxe; d++) {
                for (ptrdiff_t e = 0; e < ne; ++e) {
                    const idx_t n = elements[d][e];
                    if (n >= 0 && (size_t)n < mask.size()) {
                        mask[(size_t)n] = 1;
                    }
                }
            }
        }
    }

    static void assert_supported_ss_restrict_meshes(const Mesh &from_mesh, const Mesh &to_mesh) {
        bool has_hex  = false;
        bool has_tet  = false;
        bool has_quad = false;
        for (size_t b = 0; b < from_mesh.n_blocks(); ++b) {
            const auto bid       = static_cast<block_idx_t>(b);
            const auto from_type = from_mesh.element_type(bid);
            const auto to_type   = to_mesh.element_type(bid);
            if (!is_semistructured_type(from_type)) {
                SMESH_ERROR("Restrict: from-mesh block %zu is not semistructured (type %s)\n",
                            b,
                            type_to_string(from_type));
            }

            const auto from_family = ss_source_family(from_type);
            const auto to_family   = transfer_source_family(to_type);
            if (from_family != to_family) {
                SMESH_ERROR("Restrict: from/to block %zu family mismatch (%s vs %s)\n",
                            b,
                            type_to_string(from_type),
                            type_to_string(to_type));
            }

            if (from_family != HEX8 && from_family != QUAD4 && from_family != TET4) {
                SMESH_ERROR("Restrict: SS family %s is not implemented\n", type_to_string(from_family));
            }

            has_hex |= from_family == HEX8;
            has_tet |= from_family == TET4;
            has_quad |= from_family == QUAD4;
        }

        if (has_quad && (has_hex || has_tet)) {
            SMESH_ERROR("Restrict: mixed SS families with QUAD are not implemented\n");
        }
    }

    template <typename T>
    class Restrict<T>::Impl {
    public:
        std::shared_ptr<Mesh>     from_mesh;
        std::shared_ptr<Mesh>     to_mesh;
        ExecutionSpace            es;
        SharedBuffer<uint16_t>    element_to_node_incidence_count;
        SharedBuffer<uint16_t>    hex_incidence;
        SharedBuffer<uint16_t>    tet_incidence;
        SharedBuffer<T>           tet_to_work;
        std::vector<char>         hex_nodes_coarse;
        std::vector<char>         tet_nodes_coarse;
        bool                      mixed_hex_tet{false};
        int                       block_size;
        std::shared_ptr<Exchange> from_exchange;
        std::shared_ptr<Exchange> to_exchange;

        std::function<int(const T *const x, T *const y)> actual_op;

        static bool is_mpi_distributed(const std::shared_ptr<Mesh> &mesh) {
            return mesh && mesh->is_distributed() && mesh->comm() && mesh->comm()->size() > 1;
        }

        void init() {
            if (!from_mesh || !to_mesh) {
                SMESH_ERROR("Restrict: null mesh\n");
            }
            if (from_mesh->n_blocks() != to_mesh->n_blocks()) {
                SMESH_ERROR("Restrict: from/to n_blocks mismatch (%zu vs %zu)\n",
                            from_mesh->n_blocks(),
                            to_mesh->n_blocks());
            }

            const auto from_element = from_mesh->element_type(0);
            const auto to_element   = to_mesh->element_type(0);
            const bool from_ss      = is_semistructured_type(from_element);
            const bool to_ss        = is_semistructured_type(to_element);
            if (from_ss) {
                assert_supported_ss_restrict_meshes(*from_mesh, *to_mesh);
            }

            if (!from_ss && from_mesh->n_blocks() > 1) {
                SMESH_ERROR("Restrict: unstructured multi-block restriction is not implemented\n");
            }

            element_to_node_incidence_count = create_host_buffer<uint16_t>(from_mesh->n_nodes());
            {
                auto buff = element_to_node_incidence_count->data();
                for (size_t b = 0; b < from_mesh->n_blocks(); ++b) {
                    auto block = from_mesh->block(b);
                    // Owned elements only: aura copies would double-count after scatter-add.
                    const ptrdiff_t nelements =
                            is_mpi_distributed(from_mesh) ? block->n_elements_owned() : block->n_elements();
                    if (nelements == 0) {
                        continue;
                    }
                    const int     nxe      = block->n_nodes_per_element();
                    idx_t **const elements = block->elements()->data();
                    for (int d = 0; d < nxe; d++) {
                        for (ptrdiff_t i = 0; i < nelements; ++i) {
                            buff[elements[d][i]]++;
                        }
                    }
                }
            }

            if (is_mpi_distributed(from_mesh)) {
                from_exchange = Exchange::create_nodal(from_mesh, Exchange::ExchangeScope::GhostsAndAura);
                const ptrdiff_t nn   = from_mesh->n_nodes();
                auto            cnt  = create_host_buffer<i32>(nn);
                auto            buff = element_to_node_incidence_count->data();
                for (ptrdiff_t i = 0; i < nn; ++i) {
                    cnt->data()[i] = static_cast<i32>(buff[i]);
                }
                if (from_exchange->scatter_add(cnt->data()) != SMESH_SUCCESS ||
                    from_exchange->gather(cnt->data()) != SMESH_SUCCESS) {
                    SMESH_ERROR("Restrict: incidence scatter-add/gather failed\n");
                }
                for (ptrdiff_t i = 0; i < nn; ++i) {
                    const i32 v = cnt->data()[i];
                    buff[i]     = static_cast<uint16_t>(std::min<i32>(v, 65535));
                }
            }
            if (is_mpi_distributed(to_mesh)) {
                to_exchange = Exchange::create_nodal(to_mesh, Exchange::ExchangeScope::GhostsAndAura);
            }

            mixed_hex_tet = from_ss && ss_has_hex_and_tet(*from_mesh);
            if (mixed_hex_tet) {
                const bool mpi_from = is_mpi_distributed(from_mesh);
                auto       fill_family_incidence = [&](const ElemType family) {
                    auto buf  = create_host_buffer<uint16_t>(from_mesh->n_nodes());
                    auto buff = buf->data();
                    for (size_t b = 0; b < from_mesh->n_blocks(); ++b) {
                        if (ss_block_family(*from_mesh, b) != family) {
                            continue;
                        }
                        auto            block = from_mesh->block(b);
                        const ptrdiff_t ne    = mpi_from ? block->n_elements_owned() : block->n_elements();
                        add_block_incidence(*block, ne, buff);
                    }
                    if (mpi_from && from_exchange) {
                        const ptrdiff_t nn  = from_mesh->n_nodes();
                        auto            cnt = create_host_buffer<i32>(nn);
                        for (ptrdiff_t i = 0; i < nn; ++i) {
                            cnt->data()[i] = static_cast<i32>(buff[i]);
                        }
                        if (from_exchange->scatter_add(cnt->data()) != SMESH_SUCCESS ||
                            from_exchange->gather(cnt->data()) != SMESH_SUCCESS) {
                            SMESH_ERROR("Restrict: mixed-family incidence scatter-add/gather failed\n");
                        }
                        for (ptrdiff_t i = 0; i < nn; ++i) {
                            const i32 v = cnt->data()[i];
                            buff[i]     = static_cast<uint16_t>(std::min<i32>(v, 65535));
                        }
                    }
                    return buf;
                };
                hex_incidence = fill_family_incidence(HEX8);
                tet_incidence = fill_family_incidence(TET4);
                mark_family_nodes(*to_mesh, HEX8, hex_nodes_coarse);
                mark_family_nodes(*to_mesh, TET4, tet_nodes_coarse);
                tet_to_work = create_host_buffer<T>((ptrdiff_t)to_mesh->n_nodes() * (ptrdiff_t)block_size);
            }

#ifdef SMESH_ENABLE_CUDA
            if (EXECUTION_SPACE_DEVICE == es) {
                auto dbuff = to_device(element_to_node_incidence_count);

                if (from_ss) {
                    if (!ss_all_hex_family(*from_mesh)) {
                        SMESH_ERROR("Restrict: DEVICE SS restriction is implemented for HEX-family only\n");
                    }
                    for (size_t b = 0; b < from_mesh->n_blocks(); ++b) {
                        if (from_mesh->block(b)->n_elements() == 0) {
                            continue;
                        }
                        auto from_els = from_mesh->block(b)->device_elements_SoA();
                        if (from_els->mem_space() != MEMORY_SPACE_DEVICE) {
                            SMESH_ERROR("Elements are not on the device");
                            return;
                        }
                        if (to_ss && to_mesh->block(b)->n_elements() > 0) {
                            auto to_els = to_mesh->block(b)->device_elements_SoA();
                            if (to_els->mem_space() != MEMORY_SPACE_DEVICE) {
                                SMESH_ERROR("To elements are not on the device");
                                return;
                            }
                        }
                    }

                    if (to_ss) {
                        actual_op = [=](const real_t *const from, real_t *const to) -> int {
                            SMESH_TRACE_SCOPE("cu_sshex8_restrict");
                            const int from_level = semistructured_level(*from_mesh);
                            const int to_level   = semistructured_level(*to_mesh);
                            int       err        = SMESH_SUCCESS;
                            for (size_t b = 0; b < from_mesh->n_blocks(); ++b) {
                                auto            from_b = from_mesh->block(b);
                                auto            to_b   = to_mesh->block(b);
                                const ptrdiff_t ne     = from_b->n_elements();
                                if (ne == 0) {
                                    continue;
                                }
                                err |= smesh::cu_sshex8_restrict(ne,
                                                                 from_level,
                                                                 1,
                                                                 from_b->device_elements_SoA()->data(),
                                                                 dbuff->data(),
                                                                 to_level,
                                                                 1,
                                                                 to_b->device_elements_SoA()->data(),
                                                                 block_size,
                                                                 SMESH_DEFAULT,
                                                                 1,
                                                                 from,
                                                                 SMESH_DEFAULT,
                                                                 1,
                                                                 to,
                                                                 SMESH_DEFAULT_STREAM);
                            }
                            return err;
                        };
                        return;
                    }

                    actual_op = [=](const real_t *const from, real_t *const to) -> int {
                        SMESH_TRACE_SCOPE("cu_sshex8_hierarchical_restriction");
                        const int from_level = semistructured_level(*from_mesh);
                        int       err        = SMESH_SUCCESS;
                        for (size_t b = 0; b < from_mesh->n_blocks(); ++b) {
                            auto            from_b = from_mesh->block(b);
                            const ptrdiff_t ne     = from_b->n_elements();
                            if (ne == 0) {
                                continue;
                            }
                            err |= cu_sshex8_hierarchical_restriction(from_level,
                                                                      ne,
                                                                      from_b->device_elements_SoA()->data(),
                                                                      dbuff->data(),
                                                                      block_size,
                                                                      SMESH_DEFAULT,
                                                                      1,
                                                                      from,
                                                                      SMESH_DEFAULT,
                                                                      1,
                                                                      to,
                                                                      SMESH_DEFAULT_STREAM);
                        }
                        return err;
                    };
                    return;
                }

                auto elements = from_mesh->block(0)->device_elements_SoA();
                if (elements->mem_space() != MEMORY_SPACE_DEVICE) {
                    SMESH_ERROR("Elements are not on the device");
                    return;
                }

                actual_op = [=](const real_t *const from, real_t *const to) -> int {
                    SMESH_TRACE_SCOPE("cu_macrotet4_to_tet4_restriction_element_based");

                    return cu_macrotet4_to_tet4_restriction_element_based(from_mesh->n_elements(),
                                                                          elements->data(),
                                                                          dbuff->data(),
                                                                          block_size,
                                                                          SMESH_DEFAULT,
                                                                          1,
                                                                          from,
                                                                          SMESH_DEFAULT,
                                                                          1,
                                                                          to,
                                                                          SMESH_DEFAULT_STREAM);
                };
                return;
            } else
#endif
            {
                if (from_ss) {
                    if (!to_ss) {
                        if (is_mpi_distributed(from_mesh) || is_mpi_distributed(to_mesh)) {
                            actual_op = [=](const T *const from, T *const to) -> int {
                                SMESH_TRACE_SCOPE("mpi_ss_hierarchical_restriction");
                                const int from_level = semistructured_level(*from_mesh);
                                auto      count      = element_to_node_incidence_count->data();
                                int       err        = SMESH_SUCCESS;
                                for (size_t b = 0; b < from_mesh->n_blocks(); ++b) {
                                    auto            from_b = from_mesh->block(b);
                                    const ptrdiff_t ne     = from_b->n_elements_owned();
                                    if (ne == 0) {
                                        continue;
                                    }
                                    const auto fam = ss_block_family(*from_mesh, b);
                                    if (fam == HEX8) {
                                        auto   to_b = to_mesh->block(b);
                                        idx_t *to_sshex[8];
                                        hex8_elements_as_sshex8_level1(to_b->elements()->data(), to_sshex);
                                        err |= sshex8_restrict(ne,
                                                               from_level,
                                                               1,
                                                               from_b->elements()->data(),
                                                               count,
                                                               1,
                                                               1,
                                                               to_sshex,
                                                               block_size,
                                                               from,
                                                               to);
                                    } else if (fam == QUAD4) {
                                        err |= ssquad4_hierarchical_restriction(
                                                from_level, ne, from_b->elements()->data(), count, block_size, from, to);
                                    } else {
                                        err |= sstet4_hierarchical_restriction(
                                                from_level, ne, from_b->elements()->data(), count, block_size, from, to);
                                    }
                                }
                                return err;
                            };
                            return;
                        }

                        actual_op = [=](const T *const from, T *const to) -> int {
                            SMESH_TRACE_SCOPE("ss_hierarchical_restriction");
                            const int from_level = semistructured_level(*from_mesh);
                            auto      count      = element_to_node_incidence_count->data();
                            int       err        = SMESH_SUCCESS;
                            for (size_t b = 0; b < from_mesh->n_blocks(); ++b) {
                                auto            from_b = from_mesh->block(b);
                                const ptrdiff_t ne     = from_b->n_elements();
                                if (ne == 0) {
                                    continue;
                                }
                                const auto fam = ss_block_family(*from_mesh, b);
                                if (fam == HEX8) {
                                    err |= sshex8_hierarchical_restriction(
                                            from_level, ne, from_b->elements()->data(), count, block_size, from, to);
                                } else if (fam == QUAD4) {
                                    err |= ssquad4_hierarchical_restriction(
                                            from_level, ne, from_b->elements()->data(), count, block_size, from, to);
                                } else {
                                    err |= sstet4_hierarchical_restriction(
                                            from_level, ne, from_b->elements()->data(), count, block_size, from, to);
                                }
                            }
                            return err;
                        };
                        return;
                    }

                    actual_op = [=](const T *const from, T *const to) -> int {
                        SMESH_TRACE_SCOPE("ss_restrict");
                        const int from_level = semistructured_level(*from_mesh);
                        const int to_level   = semistructured_level(*to_mesh);
                        const int bs         = block_size;

                        auto restrict_family = [&](const ElemType family, const uint16_t *const count, T *const out) {
                            int err = SMESH_SUCCESS;
                            for (size_t b = 0; b < from_mesh->n_blocks(); ++b) {
                                if (ss_block_family(*from_mesh, b) != family) {
                                    continue;
                                }
                                auto            from_b = from_mesh->block(b);
                                auto            to_b   = to_mesh->block(b);
                                const ptrdiff_t ne     = is_mpi_distributed(from_mesh) ? from_b->n_elements_owned()
                                                                                       : from_b->n_elements();
                                if (ne == 0) {
                                    continue;
                                }
                                if (family == HEX8) {
                                    err |= sshex8_restrict(ne,
                                                           from_level,
                                                           1,
                                                           from_b->elements()->data(),
                                                           count,
                                                           to_level,
                                                           1,
                                                           to_b->elements()->data(),
                                                           bs,
                                                           from,
                                                           out);
                                } else {
                                    err |= sstet4_restrict(ne,
                                                           from_level,
                                                           1,
                                                           from_b->elements()->data(),
                                                           count,
                                                           to_level,
                                                           1,
                                                           to_b->elements()->data(),
                                                           bs,
                                                           from,
                                                           out);
                                }
                            }
                            return err;
                        };

                        if (mixed_hex_tet) {
                            auto            tet_out = tet_to_work->data();
                            const ptrdiff_t nn      = to_mesh->n_nodes();
                            std::memset(tet_out, 0, (size_t)nn * (size_t)bs * sizeof(T));
                            int err = restrict_family(HEX8, hex_incidence->data(), to);
                            err |= restrict_family(TET4, tet_incidence->data(), tet_out);
                            for (ptrdiff_t i = 0; i < nn; ++i) {
                                const bool hx = hex_nodes_coarse[(size_t)i] != 0;
                                const bool te = tet_nodes_coarse[(size_t)i] != 0;
                                for (int d = 0; d < bs; ++d) {
                                    const ptrdiff_t idx = i * (ptrdiff_t)bs + d;
                                    if (hx && te) {
                                        to[idx] = (to[idx] + tet_out[idx]) / T(2);
                                    } else if (te) {
                                        to[idx] = tet_out[idx];
                                    }
                                }
                            }
                            return err;
                        }

                        auto count = element_to_node_incidence_count->data();
                        int  err   = SMESH_SUCCESS;
                        for (size_t b = 0; b < from_mesh->n_blocks(); ++b) {
                            auto            from_b = from_mesh->block(b);
                            auto            to_b   = to_mesh->block(b);
                            const ptrdiff_t ne     = is_mpi_distributed(from_mesh) ? from_b->n_elements_owned()
                                                                                   : from_b->n_elements();
                            if (ne == 0) {
                                continue;
                            }
                            const auto fam = ss_block_family(*from_mesh, b);
                            if (fam == HEX8) {
                                err |= sshex8_restrict(ne,
                                                       from_level,
                                                       1,
                                                       from_b->elements()->data(),
                                                       count,
                                                       to_level,
                                                       1,
                                                       to_b->elements()->data(),
                                                       bs,
                                                       from,
                                                       to);
                            } else if (fam == QUAD4) {
                                err |= ssquad4_restrict(ne,
                                                        from_level,
                                                        1,
                                                        from_b->elements()->data(),
                                                        count,
                                                        to_level,
                                                        1,
                                                        to_b->elements()->data(),
                                                        bs,
                                                        from,
                                                        to);
                            } else {
                                err |= sstet4_restrict(ne,
                                                       from_level,
                                                       1,
                                                       from_b->elements()->data(),
                                                       count,
                                                       to_level,
                                                       1,
                                                       to_b->elements()->data(),
                                                       bs,
                                                       from,
                                                       to);
                            }
                        }
                        return err;
                    };
                    return;
                }

                actual_op = [=](const T *const from, T *const to) -> int {
                    SMESH_TRACE_SCOPE("hierarchical_restriction_with_counting");

                    return hierarchical_restriction(from_element,
                                                    to_element,
                                                    from_mesh->n_elements(),
                                                    from_mesh->elements(0)->data(),
                                                    element_to_node_incidence_count->data(),
                                                    block_size,
                                                    from,
                                                    to);
                };
            }
        }

        Impl(const std::shared_ptr<Mesh> &from, const std::shared_ptr<Mesh> &to, const ExecutionSpace es, const int block_size)
            : from_mesh(from), to_mesh(to), es(es), block_size(block_size) {
            init();
        }

        int apply(const T *const x, T *const y) {
            if ((from_exchange || to_exchange) && es == EXECUTION_SPACE_DEVICE) {
                SMESH_ERROR("Restrict: distributed DEVICE apply is not implemented\n");
                return SMESH_FAILURE;
            }

            T *const x_mut = const_cast<T *>(x);
            if (from_exchange) {
                if (from_exchange->gather(x_mut, block_size) != SMESH_SUCCESS) {
                    return SMESH_FAILURE;
                }
            }

            const int err = actual_op(x, y);
            if (err != SMESH_SUCCESS) {
                return err;
            }

            if (to_exchange) {
                if (to_exchange->scatter_add(y, block_size) != SMESH_SUCCESS) {
                    return SMESH_FAILURE;
                }
                if (to_exchange->gather(y, block_size) != SMESH_SUCCESS) {
                    return SMESH_FAILURE;
                }
            }
            return SMESH_SUCCESS;
        }
    };

    template <typename T>
    Restrict<T>::Restrict(const std::shared_ptr<Mesh> &from,
                          const std::shared_ptr<Mesh> &to,
                          const ExecutionSpace         es,
                          const int                    block_size)
        : impl_(std::make_unique<Impl>(from, to, es, block_size)) {}

    template <typename T>
    std::shared_ptr<Restrict<T>> Restrict<T>::create(const std::shared_ptr<Mesh> &from,
                                                     const std::shared_ptr<Mesh> &to,
                                                     const ExecutionSpace         es,
                                                     const int                    block_size) {
        return std::make_shared<Restrict<T>>(from, to, es, block_size);
    }

    template <typename T>
    Restrict<T>::~Restrict() = default;

    template <typename T>
    int Restrict<T>::apply(const T *const x, T *const y) {
        return impl_->apply(x, y);
    }

    template <typename T>
    ptrdiff_t Restrict<T>::rows() const {
        return impl_->from_mesh->n_nodes() * impl_->block_size;
    }

    template <typename T>
    ptrdiff_t Restrict<T>::cols() const {
        return impl_->to_mesh->n_nodes() * impl_->block_size;
    }

    template <typename T>
    ExecutionSpace Restrict<T>::execution_space() const {
        return impl_->es;
    }

    template <typename T>
    const SharedBuffer<uint16_t> &Restrict<T>::element_to_node_incidence_count() const {
        return impl_->element_to_node_incidence_count;
    }

    template <typename T>
    class SurfaceRestrict<T>::Impl {
    public:
        int                    from_level;
        ElemType               from_elem_type;
        ptrdiff_t              from_n_nodes;
        SharedBuffer<idx_t *>  from_sides;
        SharedBuffer<uint16_t> from_count;

        int                   to_level;
        ElemType              to_elem_type;
        ptrdiff_t             to_n_nodes;
        SharedBuffer<idx_t *> to_sides;

        ExecutionSpace es;
        int            block_size;

        std::function<int(const T *const x, T *const y)> actual_op;

        Impl(const int                     from_level,
             const ElemType                from_elem_type,
             const ptrdiff_t               from_n_nodes,
             const SharedBuffer<idx_t *>  &from_sides,
             const SharedBuffer<uint16_t> &from_count,
             const int                     to_level,
             const ElemType                to_elem_type,
             const ptrdiff_t               to_n_nodes,
             const SharedBuffer<idx_t *>  &to_sides,
             const ExecutionSpace          es,
             const int                     block_size)
            : from_level(from_level),
              from_elem_type(from_elem_type),
              from_n_nodes(from_n_nodes),
              from_sides(from_sides),
              from_count(from_count),
              to_level(to_level),
              to_elem_type(to_elem_type),
              to_n_nodes(to_n_nodes),
              to_sides(to_sides),
              es(es),
              block_size(block_size) {}

        int apply(const T *const x, T *const y) {
#ifdef SMESH_ENABLE_CUDA
            if (es == EXECUTION_SPACE_DEVICE) {
                cu_ssquad4_restrict(from_sides->extent(1),
                                    from_level,
                                    1,
                                    from_sides->data(),
                                    from_count->data(),
                                    to_level,
                                    1,
                                    to_sides->data(),
                                    block_size,
                                    SMESH_DEFAULT,
                                    1,
                                    x,
                                    SMESH_DEFAULT,
                                    1,
                                    y,
                                    SMESH_DEFAULT_STREAM);
                return SMESH_SUCCESS;
            }
#endif
            ssquad4_restrict(from_sides->extent(1),
                             from_level,
                             1,
                             from_sides->data(),
                             from_count->data(),
                             to_level,
                             1,
                             to_sides->data(),
                             block_size,
                             x,
                             y);

            return SMESH_SUCCESS;
        }
    };

    template <typename T>
    SurfaceRestrict<T>::SurfaceRestrict(const int                     from_level,
                                        const ElemType                from_elem_type,
                                        const ptrdiff_t               from_n_nodes,
                                        const SharedBuffer<idx_t *>  &from_sides,
                                        const SharedBuffer<uint16_t> &from_count,
                                        const int                     to_level,
                                        const ElemType                to_elem_type,
                                        const ptrdiff_t               to_n_nodes,
                                        const SharedBuffer<idx_t *>  &to_sides,
                                        const ExecutionSpace          es,
                                        const int                     block_size)
        : impl_(std::make_unique<Impl>(from_level,
                                       from_elem_type,
                                       from_n_nodes,
                                       from_sides,
                                       from_count,
                                       to_level,
                                       to_elem_type,
                                       to_n_nodes,
                                       to_sides,
                                       es,
                                       block_size)) {}

    template <typename T>
    std::shared_ptr<SurfaceRestrict<T>> SurfaceRestrict<T>::create(const int                     from_level,
                                                                   const ElemType                from_elem_type,
                                                                   const ptrdiff_t               from_n_nodes,
                                                                   const SharedBuffer<idx_t *>  &from_sides,
                                                                   const SharedBuffer<uint16_t> &from_count,
                                                                   const int                     to_level,
                                                                   const ElemType                to_elem_type,
                                                                   const ptrdiff_t               to_n_nodes,
                                                                   const SharedBuffer<idx_t *>  &to_sides,
                                                                   const ExecutionSpace          es,
                                                                   const int                     block_size) {
        return std::make_shared<SurfaceRestrict<T>>(from_level,
                                                    from_elem_type,
                                                    from_n_nodes,
                                                    from_sides,
                                                    from_count,
                                                    to_level,
                                                    to_elem_type,
                                                    to_n_nodes,
                                                    to_sides,
                                                    es,
                                                    block_size);
    }

    template <typename T>
    SurfaceRestrict<T>::~SurfaceRestrict() = default;

    template <typename T>
    int SurfaceRestrict<T>::apply(const T *const x, T *const y) {
        return impl_->apply(x, y);
    }

    template <typename T>
    ptrdiff_t SurfaceRestrict<T>::rows() const {
        return impl_->to_n_nodes * impl_->block_size;
    }

    template <typename T>
    ptrdiff_t SurfaceRestrict<T>::cols() const {
        return impl_->from_n_nodes * impl_->block_size;
    }

    template <typename T>
    ExecutionSpace SurfaceRestrict<T>::execution_space() const {
        return impl_->es;
    }

    template class Restrict<real_t>;
    template class SurfaceRestrict<real_t>;

}  // namespace smesh

