#include "smesh_restrict.hpp"
#include "smesh_buffer.hpp"
#include "smesh_mesh.hpp"
#include "smesh_restriction.hpp"
#include "smesh_semistructured.hpp"
#include "smesh_sshex8.hpp"
#include "smesh_ssquad4.hpp"
#include "smesh_sshex8_restriction.hpp"
#include "smesh_ssquad4_restriction.hpp"
#include "smesh_tracer.hpp"

#ifdef SMESH_ENABLE_CUDA
#include "cu_sshex8_interpolate.hpp"
#include "cu_ssquad4_interpolate.hpp"
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

        std::function<int(const T* const x, T* const y)> actual_op;

        void init() {
            auto from_element = from_mesh->element_type(0);
            auto to_element   = to_mesh->element_type(0);

            ptrdiff_t nnodes   = 0;
            idx_t**   elements = nullptr;
            int       nxe;
            if (is_semistructured_type(from_element)) {
                nxe      = sshex8_nxe(semistructured_level(*from_mesh));
                elements = from_mesh->elements(0)->data();
                nnodes   = from_mesh->n_nodes();
            } else {
                nxe      = elem_num_nodes(from_element);
                elements = from_mesh->elements(0)->data();
                nnodes   = from_mesh->n_nodes();
            }

            element_to_node_incidence_count = create_host_buffer<uint16_t>(nnodes);
            {
                auto buff = element_to_node_incidence_count->data();

                // #pragma omp parallel for // BAD performance with parallel for
                const ptrdiff_t nelements = from_mesh->n_elements();
                for (int d = 0; d < nxe; d++) {
                    for (ptrdiff_t i = 0; i < nelements; ++i) {
                        // #pragma omp atomic update
                        buff[elements[d][i]]++;
                    }
                }
            }

#ifdef SMESH_ENABLE_CUDA
            if (EXECUTION_SPACE_DEVICE == es) {
                auto dbuff = to_device(element_to_node_incidence_count);

                auto elements = from_mesh->device_elements();
                if (!elements) {
                    elements = create_device_elements(from_mesh, from_mesh->element_type());
                    from_mesh->set_device_elements(elements);
                }

                if (is_semistructured_type(from_element)) {
                    if (is_semistructured_type(to_element)) {
                        // FIXME make sure to reuse fine level elements and strides
                        auto to_elements = to_mesh->device_elements();
                        if (!to_elements) {
                            to_elements = create_device_elements(to_mesh, to_mesh->element_type());
                            to_mesh->set_device_elements(to_elements);
                        }

                        actual_op = [=](const real_t* const from, real_t* const to) -> int {
                            SMESH_TRACE_SCOPE("cu_sshex8_restrict");
                            return cu_sshex8_restrict(from_mesh->n_elements(),
                                                      semistructured_level(*from_mesh),
                                                      1,
                                                      elements->data(),
                                                      dbuff->data(),
                                                      semistructured_level(*to_mesh),
                                                      1,
                                                      to_elements->data(),
                                                      block_size,
                                                      SMESH_REAL_DEFAULT,
                                                      1,
                                                      from,
                                                      SMESH_REAL_DEFAULT,
                                                      1,
                                                      to,
                                                      SMESH_DEFAULT_STREAM);
                        };
                        return;
                    } else {
                        actual_op = [=](const real_t* const from, real_t* const to) -> int {
                            SMESH_TRACE_SCOPE("cu_sshex8_hierarchical_restriction");

                            return cu_sshex8_hierarchical_restriction(semistructured_level(*from_mesh),
                                                                      from_mesh->n_elements(),
                                                                      elements->data(),
                                                                      dbuff->data(),
                                                                      block_size,
                                                                      SMESH_REAL_DEFAULT,
                                                                      1,
                                                                      from,
                                                                      SMESH_REAL_DEFAULT,
                                                                      1,
                                                                      to,
                                                                      SMESH_DEFAULT_STREAM);
                        };
                        return;
                    }
                } else {
                    actual_op = [=](const real_t* const from, real_t* const to) -> int {
                        SMESH_TRACE_SCOPE("cu_macrotet4_to_tet4_restriction_element_based");

                        return cu_macrotet4_to_tet4_restriction_element_based(from_mesh->n_elements(),
                                                                              elements->data(),
                                                                              dbuff->data(),
                                                                              block_size,
                                                                              SMESH_REAL_DEFAULT,
                                                                              1,
                                                                              from,
                                                                              SMESH_REAL_DEFAULT,
                                                                              1,
                                                                              to,
                                                                              SMESH_DEFAULT_STREAM);
                    };
                    return;
                }
            } else
#endif
            {
                if (is_semistructured_type(from_element)) {
                    if (!is_semistructured_type(to_element)) {
                        actual_op = [=](const real_t* const from, real_t* const to) -> int {
                            SMESH_TRACE_SCOPE("sshex8_hierarchical_restriction");

                            return sshex8_hierarchical_restriction(semistructured_level(*from_mesh),
                                                                   from_mesh->n_elements(),
                                                                   from_mesh->elements(0)->data(),
                                                                   element_to_node_incidence_count->data(),
                                                                   block_size,
                                                                   from,
                                                                   to);
                        };
                        return;
                    } else {
                        actual_op = [=](const real_t* const from, real_t* const to) -> int {
                            SMESH_TRACE_SCOPE("sshex8_restrict");

                            return sshex8_restrict(from_mesh->n_elements(),           // nelements,
                                                   semistructured_level(*from_mesh),  // from_level
                                                   1,                                 // from_level_stride
                                                   from_mesh->elements(0)->data(),    // from_elements
                                                   element_to_node_incidence_count->data(),
                                                   semistructured_level(*to_mesh),  // to_level
                                                   1,                               // to_level_stride
                                                   to_mesh->elements(0)->data(),    // to_elements
                                                   block_size,                      // vec_size
                                                   from,
                                                   to);
                        };
                        return;
                    }
                } else {
                    actual_op = [=](const real_t* const from, real_t* const to) -> int {
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
                    return;
                }
            }
        }

        Impl(const std::shared_ptr<Mesh>& from, const std::shared_ptr<Mesh>& to, const ExecutionSpace es, const int block_size)
            : from_mesh(from), to_mesh(to), es(es), block_size(block_size) {
            init();
        }

        int apply(const T* const x, T* const y) { return actual_op(x, y); }
    };

    template <typename T>
    Restrict<T>::Restrict(const std::shared_ptr<Mesh>& from,
                          const std::shared_ptr<Mesh>& to,
                          const ExecutionSpace         es,
                          const int                    block_size)
        : impl_(std::make_unique<Impl>(from, to, es, block_size)) {}

    template <typename T>
    std::shared_ptr<Restrict<T>> Restrict<T>::create(const std::shared_ptr<Mesh>& from,
                                                     const std::shared_ptr<Mesh>& to,
                                                     const ExecutionSpace         es,
                                                     const int                    block_size) {
        return std::make_shared<Restrict<T>>(from, to, es, block_size);
    }

    template <typename T>
    Restrict<T>::~Restrict() = default;

    template <typename T>
    int Restrict<T>::apply(const T* const x, T* const y) {
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
    const SharedBuffer<uint16_t>& Restrict<T>::element_to_node_incidence_count() const {
        return impl_->element_to_node_incidence_count;
    }

    template <typename T>
    class SurfaceRestrict<T>::Impl {
    public:
        int                    from_level;
        smesh::ElemType        from_elem_type;
        ptrdiff_t              from_n_nodes;
        SharedBuffer<idx_t*>   from_sides;
        SharedBuffer<uint16_t> from_count;

        int                  to_level;
        smesh::ElemType      to_elem_type;
        ptrdiff_t            to_n_nodes;
        SharedBuffer<idx_t*> to_sides;

        ExecutionSpace es;
        int            block_size;

        std::function<int(const T* const x, T* const y)> actual_op;

        Impl(const int                     from_level,
             const smesh::ElemType         from_elem_type,
             const ptrdiff_t               from_n_nodes,
             const SharedBuffer<idx_t*>&   from_sides,
             const SharedBuffer<uint16_t>& from_count,
             const int                     to_level,
             const smesh::ElemType         to_elem_type,
             const ptrdiff_t               to_n_nodes,
             const SharedBuffer<idx_t*>&   to_sides,
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

        int apply(const T* const x, T* const y) {
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
                                    SMESH_REAL_DEFAULT,
                                    1,
                                    x,
                                    SMESH_REAL_DEFAULT,
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
                                        const smesh::ElemType         from_elem_type,
                                        const ptrdiff_t               from_n_nodes,
                                        const SharedBuffer<idx_t*>&   from_sides,
                                        const SharedBuffer<uint16_t>& from_count,
                                        const int                     to_level,
                                        const smesh::ElemType         to_elem_type,
                                        const ptrdiff_t               to_n_nodes,
                                        const SharedBuffer<idx_t*>&   to_sides,
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
                                                                   const smesh::ElemType         from_elem_type,
                                                                   const ptrdiff_t               from_n_nodes,
                                                                   const SharedBuffer<idx_t*>&   from_sides,
                                                                   const SharedBuffer<uint16_t>& from_count,
                                                                   const int                     to_level,
                                                                   const smesh::ElemType         to_elem_type,
                                                                   const ptrdiff_t               to_n_nodes,
                                                                   const SharedBuffer<idx_t*>&   to_sides,
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
    int SurfaceRestrict<T>::apply(const T* const x, T* const y) {
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
