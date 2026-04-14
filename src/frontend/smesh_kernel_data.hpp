#ifndef SMESH_KERNEL_DATA_HPP
#define SMESH_KERNEL_DATA_HPP

#include "smesh_buffer.hpp"
#include "smesh_forward_declarations.hpp"

#include <memory>

namespace smesh {
    class Points {
    public:
        static std::shared_ptr<Points> create_AoS(const SharedBuffer<geom_t *> &points, const MemorySpace space);
        static std::shared_ptr<Points> create_SoA(const SharedBuffer<geom_t *> &points, const MemorySpace space);

        int init_AoS(const SharedBuffer<geom_t *> &points, const MemorySpace space);
        int init_SoA(const SharedBuffer<geom_t *> &points, const MemorySpace space);

        inline SharedBuffer<geom_t *> points_SoA() const {
            SMESH_ASSERT(points_SoA_);
            return points_SoA_;
        }

        inline SharedBuffer<geom_t> points_AoS() const {
            SMESH_ASSERT(points_AoS_);
            return points_AoS_;
        }

        bool has_SoA() const { return points_SoA_ != nullptr; }
        bool has_AoS() const { return points_AoS_ != nullptr; }


    private:
        SharedBuffer<geom_t *> points_SoA_;
        SharedBuffer<geom_t>   points_AoS_;
    };

    class Elements {
    public:
        static std::shared_ptr<Elements> create_AoS(const SharedBuffer<idx_t *> &elements, const MemorySpace space);
        static std::shared_ptr<Elements> create_SoA(const SharedBuffer<idx_t *> &elements, const MemorySpace space);

        int init_AoS(const SharedBuffer<idx_t *> &elements, const MemorySpace space);
        int init_SoA(const SharedBuffer<idx_t *> &elements, const MemorySpace space);

        inline SharedBuffer<idx_t *> elements_SoA() const {
            SMESH_ASSERT(elements_SoA_);
            return elements_SoA_;
        }

        inline SharedBuffer<idx_t> elements_AoS() const {
            SMESH_ASSERT(elements_AoS_);
            return elements_AoS_;
        }

        bool has_SoA() const { return elements_SoA_ != nullptr; }
        bool has_AoS() const { return elements_AoS_ != nullptr; }

    private:
        SharedBuffer<idx_t *> elements_SoA_;
        SharedBuffer<idx_t>   elements_AoS_;
    };

    class Jacobian {
    public:
        static std::shared_ptr<Jacobian> create_AoS(const std::shared_ptr<Mesh> &mesh,
                                                    const MemorySpace            space,
                                                    const block_idx_t            block_id);

        static std::shared_ptr<Jacobian> create_SoA(const std::shared_ptr<Mesh> &mesh,
                                                    const MemorySpace            space,
                                                    const block_idx_t            block_id);

        inline SharedBuffer<jacobian_t *> jacobians_SoA() const {
            SMESH_ASSERT(jacobians_SoA_);
            return jacobians_SoA_;
        }

        inline SharedBuffer<jacobian_t> jacobians_AoS() const {
            SMESH_ASSERT(jacobians_AoS_);
            return jacobians_AoS_;
        }

    private:
        SharedBuffer<jacobian_t *> jacobians_SoA_;
        SharedBuffer<jacobian_t>   jacobians_AoS_;
    };

    class JacobianAdjugateAndDeterminant {
    public:
        static std::shared_ptr<JacobianAdjugateAndDeterminant> create_AoS(const std::shared_ptr<Mesh> &mesh,
                                                                          const MemorySpace            space,
                                                                          const block_idx_t            block_id);

        static std::shared_ptr<JacobianAdjugateAndDeterminant> create_SoA(const std::shared_ptr<Mesh> &mesh,
                                                                          const MemorySpace            space,
                                                                          const block_idx_t            block_id);

        inline SharedBuffer<jacobian_t *> jacobian_adjugate_SoA() const {
            SMESH_ASSERT(jacobian_adjugate_SoA_);
            return jacobian_adjugate_SoA_;
        }

        inline SharedBuffer<jacobian_t> jacobian_adjugate_AoS() const {
            SMESH_ASSERT(jacobian_adjugate_AoS_);
            return jacobian_adjugate_AoS_;
        }

        inline SharedBuffer<geom_t> jacobian_determinant() const {
            SMESH_ASSERT(jacobian_determinant_);
            return jacobian_determinant_;
        }

    private:
        SharedBuffer<jacobian_t *> jacobian_adjugate_SoA_;
        SharedBuffer<jacobian_t>   jacobian_adjugate_AoS_;
        SharedBuffer<geom_t>       jacobian_determinant_;
    };

    class FFF {
    public:
        static std::shared_ptr<FFF> create_AoS(const std::shared_ptr<Mesh> &mesh,
                                               const MemorySpace            space,
                                               const block_idx_t            block_id);

        static std::shared_ptr<FFF> create_SoA(const std::shared_ptr<Mesh> &mesh,
                                               const MemorySpace            space,
                                               const block_idx_t            block_id);

        inline SharedBuffer<jacobian_t *> fff_SoA() const {
            SMESH_ASSERT(fff_SoA_);
            return fff_SoA_;
        }

        inline SharedBuffer<jacobian_t> fff_AoS() const {
            SMESH_ASSERT(fff_AoS_);
            return fff_AoS_;
        }

    private:
        SharedBuffer<jacobian_t *> fff_SoA_;
        SharedBuffer<jacobian_t>   fff_AoS_;
    };
}  // namespace smesh

#endif  // SMESH_KERNEL_DATA_HPP
