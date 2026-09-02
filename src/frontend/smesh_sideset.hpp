#ifndef SMESH_SIDESET_HPP
#define SMESH_SIDESET_HPP

#include "smesh_base.hpp"
#include "smesh_buffer.hpp"
#include "smesh_forward_declarations.hpp"

#include <cstddef>
#include <functional>
#include <memory>
#include <vector>

namespace smesh {

    /// Lightweight `(block_id, parent[], lfi[])` side set.
    /// Semistructured level changes keep the same `(parent, lfi)` valid (macro
    /// element and local face index are unchanged). Unstructured refine/promote
    /// does not remap sidesets.
    class Sideset final {
    public:
        int read(const std::shared_ptr<Communicator> &comm, const Path &path, block_idx_t block_id = 0);

        int read_and_redistibute(const std::shared_ptr<Mesh> &mesh, const Path &path, block_idx_t block_id = 0);
        int redistribute(const std::shared_ptr<Mesh> &mesh);

        SharedBuffer<element_idx_t>     parent() const;
        SharedBuffer<i16>               lfi() const;
        block_idx_t                     block_id() const;
        SharedBuffer<large_idx_t>       element_mapping() const;
        static std::shared_ptr<Sideset> create_from_file(const std::shared_ptr<Communicator> &comm,
                                                         const Path                          &path,
                                                         block_idx_t                          block_id = 0);

        ptrdiff_t                     size() const;
        std::shared_ptr<Communicator> comm() const;
        int                           write(const Path &path) const;

        Sideset(const std::shared_ptr<Communicator> &comm,
                const SharedBuffer<element_idx_t>   &parent,
                const SharedBuffer<i16>             &lfi,
                block_idx_t                          block_id        = 0,
                const SharedBuffer<large_idx_t>     &element_mapping = nullptr);
        Sideset();
        ~Sideset();

        static std::shared_ptr<Sideset> create(const std::shared_ptr<Communicator> &comm,
                                               const SharedBuffer<element_idx_t>   &parent,
                                               const SharedBuffer<i16>             &lfi,
                                               block_idx_t                          block_id        = 0,
                                               const SharedBuffer<large_idx_t>     &element_mapping = nullptr);

        static std::vector<std::shared_ptr<Sideset>> create_from_selector(
                const std::shared_ptr<Mesh>                                         &mesh,
                const std::function<bool(const geom_t, const geom_t, const geom_t)> &selector,
                const std::vector<std::string>                                      &block_names = {});

        // Preferred for performance
        static std::vector<std::shared_ptr<Sideset>> create_from_batch_selector(
                const std::shared_ptr<Mesh> &mesh,
                const std::function<
                        void(const ptrdiff_t, const geom_t *const, const geom_t *const, const geom_t *const, u8 *const selected)>
                                               &selector,
                const std::vector<std::string> &block_names = {});

        static std::vector<std::shared_ptr<Sideset>> create_from_plane(const std::shared_ptr<Mesh>    &mesh,
                                                                       const geom_t                    normal_x,
                                                                       const geom_t                    normal_y,
                                                                       const geom_t                    normal_z,
                                                                       const geom_t                    distance,
                                                                       const geom_t                    tol         = 1e-6,
                                                                       const std::vector<std::string> &block_names = {});

        /// `parent[i] <- old_to_new[parent[i]]`. `n` is the local element count.
        int remap_parents(const element_idx_t *old_to_new, ptrdiff_t n);

        /// `new_to_old[new] = old` (gather permutation from element reorder).
        int apply_element_permutation(const element_idx_t *new_to_old, ptrdiff_t n);

        void print(std::ostream &os = std::cout) const;

    private:
        class Impl;
        std::unique_ptr<Impl> impl_;
    };

    std::pair<enum ElemType, std::shared_ptr<Buffer<idx_t *>>> create_surface_from_sideset(
            const std::shared_ptr<Mesh>    &mesh,
            const std::shared_ptr<Sideset> &sideset);

    std::pair<enum ElemType, std::shared_ptr<Buffer<idx_t *>>> create_surface_from_sidesets(
            const std::shared_ptr<Mesh>                 &mesh,
            const std::vector<std::shared_ptr<Sideset>> &sideset);

    std::vector<std::pair<enum ElemType, std::shared_ptr<Buffer<idx_t *>>>> create_surfaces_from_sidesets(
            const std::shared_ptr<Mesh>                 &mesh,
            const std::vector<std::shared_ptr<Sideset>> &sidesets);

    /// Split a WEDGE/PYRAMID (or SS equivalent) sideset into homogeneous TRI and QUAD
    /// sidesets. Homogeneous or empty input is returned unchanged as a 1-element vector.
    std::vector<std::shared_ptr<Sideset>> split_mixed_arity_sideset(const std::shared_ptr<Mesh>    &mesh,
                                                                    const std::shared_ptr<Sideset> &sideset);

    SharedBuffer<idx_t> create_nodeset_from_sideset(const std::shared_ptr<Mesh> &mesh, const std::shared_ptr<Sideset> &sideset);

    SharedBuffer<idx_t> create_nodeset_from_sidesets(const std::shared_ptr<Mesh>                 &mesh,
                                                     const std::vector<std::shared_ptr<Sideset>> &sidesets);

    int remap_sidesets(const std::vector<std::shared_ptr<Sideset>> &sidesets,
                       block_idx_t                                  block_id,
                       const element_idx_t                         *old_to_new,
                       ptrdiff_t                                    n);

    /// Expand a sideset through unstructured `refine()`. SS level changes keep
    /// `(parent, lfi)` and must not use this. Promote with unchanged element count
    /// is identity. Returns nullptr (stderr diagnostic) for unsupported type pairs.
    std::shared_ptr<Sideset> map_sideset_through_refine(const std::shared_ptr<Mesh>    &coarse_mesh,
                                                        const std::shared_ptr<Sideset> &coarse_ss,
                                                        const std::shared_ptr<Mesh>    &fine_mesh);

}  // namespace smesh

#endif  // SMESH_SIDESET_HPP