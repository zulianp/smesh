#include "smesh_kernel_data.hpp"
#include "smesh_device_buffer.hpp"
#include "smesh_elem_type.hpp"
#include "smesh_fff.hpp"
#include "smesh_jacobians.hpp"
#include "smesh_mesh.hpp"

namespace smesh {
    namespace {
        template <typename T>
        SharedBuffer<T> to_memory_space(const SharedBuffer<T> &buffer, const MemorySpace space) {
            if (!buffer) {
                return nullptr;
            }

            if (space == MEMORY_SPACE_DEVICE) {
                return to_device(buffer);
            }

            return to_host(buffer);
        }

        template <typename T>
        SharedBuffer<T *> to_memory_space(const SharedBuffer<T *> &buffer, const MemorySpace space) {
            if (!buffer) {
                return nullptr;
            }

            if (space == MEMORY_SPACE_DEVICE) {
                return to_device(buffer);
            }

            return to_host(buffer);
        }

        template <typename T>
        SharedBuffer<T> soa_to_element_major_aos(const SharedBuffer<T *> &buffer) {
            return soa_to_aos(1, buffer->extent(0), buffer);
        }

        int fill_fff(const enum ElemType                                      element_type,
                     const ptrdiff_t                                          nelements,
                     const idx_t *const SMESH_RESTRICT *const SMESH_RESTRICT  elements,
                     const geom_t *const SMESH_RESTRICT *const SMESH_RESTRICT points,
                     const ptrdiff_t                                          stride,
                     jacobian_t *const SMESH_RESTRICT *const SMESH_RESTRICT   fff) {
            if (is_semistructured_type(element_type)) {
                return sshex8_macro_fff_fill(
                        proteus_hex_micro_elements_per_dim(element_type), nelements, elements, points, stride, fff);
            }

            switch (element_type) {
                case TET4:
                    return tet4_fff_fill(nelements, elements, points, stride, fff);
                case HEX8:
                    return hex8_fff_fill(nelements, elements, points, 0.5, 0.5, 0.5, stride, fff);
                default:
                    SMESH_ERROR("fill_fff: Unsupported element type: %s\n", type_to_string(element_type));
                    return SMESH_FAILURE;
            }
        }
    }  // namespace

    std::shared_ptr<Points> Points::create_AoS(const SharedBuffer<geom_t *> &points, const MemorySpace space) {
        if (!points) {
            return nullptr;
        }

        auto ret = std::make_shared<Points>();
        ret->init_AoS(points, space);
        return ret;
    }

    int Points::init_AoS(const SharedBuffer<geom_t *> &points, const MemorySpace space) {
        if (!points) {
            return SMESH_FAILURE;
        }

        auto host_points = to_host(points);
        points_AoS_      = to_memory_space(soa_to_element_major_aos(host_points), space);
        return SMESH_SUCCESS;
    }

    std::shared_ptr<Points> Points::create_SoA(const SharedBuffer<geom_t *> &points, const MemorySpace space) {
        if (!points) {
            return nullptr;
        }

        auto ret = std::make_shared<Points>();
        ret->init_SoA(points, space);
        return ret;
    }

    int Points::init_SoA(const SharedBuffer<geom_t *> &points, const MemorySpace space) {
        if (!points) {
            return SMESH_FAILURE;
        }

        points_SoA_ = to_memory_space(points, space);
        return SMESH_SUCCESS;
    }

    std::shared_ptr<Elements> Elements::create_AoS(const SharedBuffer<idx_t *> &elements, const MemorySpace space) {
        if (!elements) {
            return nullptr;
        }

        auto ret = std::make_shared<Elements>();
        ret->init_AoS(elements, space);
        return ret;
    }

    int Elements::init_AoS(const SharedBuffer<idx_t *> &elements, const MemorySpace space) {
        if (!elements) {
            return SMESH_FAILURE;
        }

        auto host_elements = to_host(elements);
        elements_AoS_      = to_memory_space(soa_to_element_major_aos(host_elements), space);
        return SMESH_SUCCESS;
    }

    std::shared_ptr<Elements> Elements::create_SoA(const SharedBuffer<idx_t *> &elements, const MemorySpace space) {
        if (!elements) {
            return nullptr;
        }

        auto ret = std::make_shared<Elements>();
        ret->init_SoA(elements, space);
        return ret;
    }

    int Elements::init_SoA(const SharedBuffer<idx_t *> &elements, const MemorySpace space) {
        if (!elements) {
            return SMESH_FAILURE;
        }

        elements_SoA_ = to_memory_space(elements, space);
        return SMESH_SUCCESS;
    }

    std::shared_ptr<Jacobian> Jacobian::create_AoS(const std::shared_ptr<Mesh> & /*mesh*/,
                                                   const MemorySpace /*space*/,
                                                   const block_idx_t /*block_id*/) {
        SMESH_ERROR("Jacobian::create_AoS is not supported yet!\n");
        return nullptr;
    }

    std::shared_ptr<Jacobian> Jacobian::create_SoA(const std::shared_ptr<Mesh> & /*mesh*/,
                                                   const MemorySpace /*space*/,
                                                   const block_idx_t /*block_id*/) {
        SMESH_ERROR("Jacobian::create_SoA is not supported yet!\n");
        return nullptr;
    }

    std::shared_ptr<JacobianAdjugateAndDeterminant> JacobianAdjugateAndDeterminant::create_AoS(const std::shared_ptr<Mesh> &mesh,
                                                                                               const MemorySpace            space,
                                                                                               const block_idx_t block_id) {
        if (!mesh) {
            SMESH_ERROR("JacobianAdjugateAndDeterminant::create_AoS: mesh is nullptr\n");
            return nullptr;
        }

        constexpr ptrdiff_t adjugate_size = 9;

        auto       ret           = std::make_shared<JacobianAdjugateAndDeterminant>();
        auto       host_elements = to_host(mesh->elements(block_id));
        auto       host_points   = to_host(mesh->points());
        const auto nelems        = mesh->n_elements(block_id);

        auto adjugate          = create_host_buffer<jacobian_t>(adjugate_size * nelems);
        auto adjugate_fake_soa = convert_host_buffer_to_fake_SoA(adjugate_size, adjugate);
        auto determinant       = create_host_buffer<geom_t>(nelems);

        if (adjugate_fill(mesh->element_type(block_id),
                          nelems,
                          host_elements->data(),
                          host_points->data(),
                          adjugate_size,
                          adjugate_fake_soa->data(),
                          determinant->data()) != SMESH_SUCCESS) {
            return nullptr;
        }

        ret->jacobian_adjugate_AoS_ = to_memory_space(adjugate, space);
        ret->jacobian_determinant_  = to_memory_space(determinant, space);
        return ret;
    }

    std::shared_ptr<JacobianAdjugateAndDeterminant> JacobianAdjugateAndDeterminant::create_SoA(const std::shared_ptr<Mesh> &mesh,
                                                                                               const MemorySpace            space,
                                                                                               const block_idx_t block_id) {
        if (!mesh) {
            SMESH_ERROR("JacobianAdjugateAndDeterminant::create_SoA: mesh is nullptr\n");
            return nullptr;
        }

        constexpr ptrdiff_t adjugate_size = 9;

        auto       ret           = std::make_shared<JacobianAdjugateAndDeterminant>();
        auto       host_elements = to_host(mesh->elements(block_id));
        auto       host_points   = to_host(mesh->points());
        const auto nelems        = mesh->n_elements(block_id);

        auto adjugate    = create_host_buffer<jacobian_t>(adjugate_size, nelems);
        auto determinant = create_host_buffer<geom_t>(nelems);

        if (adjugate_fill(mesh->element_type(block_id),
                          nelems,
                          host_elements->data(),
                          host_points->data(),
                          1,
                          adjugate->data(),
                          determinant->data()) != SMESH_SUCCESS) {
            return nullptr;
        }

        ret->jacobian_adjugate_SoA_ = to_memory_space(adjugate, space);
        ret->jacobian_determinant_  = to_memory_space(determinant, space);
        return ret;
    }

    std::shared_ptr<FFF> FFF::create_AoS(const std::shared_ptr<Mesh> &mesh, const MemorySpace space, const block_idx_t block_id) {
        if (!mesh) {
            SMESH_ERROR("FFF::create_AoS: mesh is nullptr\n");
            return nullptr;
        }

        constexpr ptrdiff_t fff_size = 6;

        auto       ret           = std::make_shared<FFF>();
        auto       host_elements = to_host(mesh->elements(block_id));
        auto       host_points   = to_host(mesh->points());
        const auto nelems        = mesh->n_elements(block_id);

        auto fff          = create_host_buffer<jacobian_t>(fff_size * nelems);
        auto fff_fake_soa = convert_host_buffer_to_fake_SoA(fff_size, fff);

        if (fill_fff(mesh->element_type(block_id),
                     nelems,
                     host_elements->data(),
                     host_points->data(),
                     fff_size,
                     fff_fake_soa->data()) != SMESH_SUCCESS) {
            return nullptr;
        }

        ret->fff_AoS_ = to_memory_space(fff, space);
        return ret;
    }

    std::shared_ptr<FFF> FFF::create_SoA(const std::shared_ptr<Mesh> &mesh, const MemorySpace space, const block_idx_t block_id) {
        if (!mesh) {
            SMESH_ERROR("FFF::create_SoA: mesh is nullptr\n");
            return nullptr;
        }

        constexpr ptrdiff_t fff_size = 6;

        auto       ret           = std::make_shared<FFF>();
        auto       host_elements = to_host(mesh->elements(block_id));
        auto       host_points   = to_host(mesh->points());
        const auto nelems        = mesh->n_elements(block_id);

        auto fff = create_host_buffer<jacobian_t>(fff_size, nelems);

        if (fill_fff(mesh->element_type(block_id), nelems, host_elements->data(), host_points->data(), 1, fff->data()) !=
            SMESH_SUCCESS) {
            return nullptr;
        }

        ret->fff_SoA_ = to_memory_space(fff, space);
        return ret;
    }

}  // namespace smesh
