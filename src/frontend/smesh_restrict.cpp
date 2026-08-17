#include "smesh_restrict.hpp"
#include "smesh_buffer.hpp"
#include "smesh_device_buffer.hpp"
#include "smesh_mesh.hpp"
#include "smesh_restriction.hpp"
#include "smesh_semistructured.hpp"
#include "smesh_sshex8.hpp"
#include "smesh_sshex8_restriction.hpp"
#include "smesh_ssquad4.hpp"
#include "smesh_ssquad4_restriction.hpp"
#include "smesh_tracer.hpp"

#ifdef SMESH_ENABLE_CUDA
#include "smesh_sshex8_restriction.cuh"
#include "smesh_ssquad4_restriction.cuh"
#include "smesh_tet4_prolongation_restriction.cuh"
#endif

namespace smesh {

    template <typename T>
    class Restrict<T>::Impl {
    public:
        std::shared_ptr<Mesh>  from_mesh;
        std::shared_ptr<Mesh>  to_mesh;
        ExecutionSpace         es;
        SharedBuffer<uint16_t> element_to_node_incidence_count;
        int                    block_size;

        std::function<int(const T *const x, T *const y)> actual_op;

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
                for (size_t b = 0; b < from_mesh->n_blocks(); ++b) {
                    const auto bid = static_cast<block_idx_t>(b);
                    if (!is_hex_ss_family(from_mesh->element_type(bid))) {
                        SMESH_ERROR(
                                "Restrict: SS restriction is HEX-family only (TET/QUAD B5.5, mixed B5.6); "
                                "from block %zu type %s\n",
                                b,
                                type_to_string(from_mesh->element_type(bid)));
                    }
                    if (!is_hex_ss_family(to_mesh->element_type(bid))) {
                        SMESH_ERROR("Restrict: to-mesh block %zu is not HEX-family (type %s)\n",
                                    b,
                                    type_to_string(to_mesh->element_type(bid)));
                    }
                }
            } else if (from_mesh->n_blocks() > 1) {
                SMESH_ERROR("Restrict: unstructured multi-block restriction is not implemented\n");
            }

            element_to_node_incidence_count = create_host_buffer<uint16_t>(from_mesh->n_nodes());
            {
                auto buff = element_to_node_incidence_count->data();
                for (size_t b = 0; b < from_mesh->n_blocks(); ++b) {
                    auto            block     = from_mesh->block(b);
                    const ptrdiff_t nelements = block->n_elements();
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

#ifdef SMESH_ENABLE_CUDA
            if (EXECUTION_SPACE_DEVICE == es) {
                auto dbuff = to_device(element_to_node_incidence_count);

                if (from_ss) {
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
                        actual_op = [=](const T *const from, T *const to) -> int {
                            SMESH_TRACE_SCOPE("sshex8_hierarchical_restriction");
                            const int from_level = semistructured_level(*from_mesh);
                            auto      count      = element_to_node_incidence_count->data();
                            int       err        = SMESH_SUCCESS;
                            for (size_t b = 0; b < from_mesh->n_blocks(); ++b) {
                                auto            from_b = from_mesh->block(b);
                                const ptrdiff_t ne     = from_b->n_elements();
                                if (ne == 0) {
                                    continue;
                                }
                                err |= sshex8_hierarchical_restriction(
                                        from_level, ne, from_b->elements()->data(), count, block_size, from, to);
                            }
                            return err;
                        };
                        return;
                    }

                    actual_op = [=](const T *const from, T *const to) -> int {
                        SMESH_TRACE_SCOPE("sshex8_restrict");
                        const int from_level = semistructured_level(*from_mesh);
                        const int to_level   = semistructured_level(*to_mesh);
                        auto      count      = element_to_node_incidence_count->data();
                        int       err        = SMESH_SUCCESS;
                        for (size_t b = 0; b < from_mesh->n_blocks(); ++b) {
                            auto            from_b = from_mesh->block(b);
                            auto            to_b   = to_mesh->block(b);
                            const ptrdiff_t ne     = from_b->n_elements();
                            if (ne == 0) {
                                continue;
                            }
                            err |= sshex8_restrict(ne,
                                                   from_level,
                                                   1,
                                                   from_b->elements()->data(),
                                                   count,
                                                   to_level,
                                                   1,
                                                   to_b->elements()->data(),
                                                   block_size,
                                                   from,
                                                   to);
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

        int apply(const T *const x, T *const y) { return actual_op(x, y); }
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

