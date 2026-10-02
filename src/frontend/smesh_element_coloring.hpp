#ifndef SMESH_ELEMENT_COLORING_HPP
#define SMESH_ELEMENT_COLORING_HPP

#include "smesh_base.hpp"

#include "smesh_buffer.hpp"
#include "smesh_forward_declarations.hpp"
#include "smesh_mesh.hpp"

#include <memory>
#include <string>
#include <vector>

namespace smesh {

    /// Colours the elements of each block so that no two elements of a colour share a node.
    ///
    /// A sweep over one colour can therefore write the global arrays with a plain `+=`: there is
    /// no conflict to resolve, so no atomics and no colour-private buffer, and the only
    /// synchronisation is one barrier between colours. The summation order a node sees is then
    /// fixed by the colouring -- one contribution per colour, in colour order -- which makes an
    /// operator built this way reproducible as long as the colouring is, and it is: the algorithm
    /// below is serial and deterministic for a fixed mesh and element order.
    ///
    /// With `modify_mesh = true` the elements are RENUMBERED into colour order, so colour `c` of
    /// block `b` is the contiguous range `[color_ptr(b)[c], color_ptr(b)[c + 1])` and a sweep
    /// needs no indirection at all. This rewrites the mesh's element arrays in place, exactly as
    /// `PackedMesh::create(..., modify_mesh = true)` rewrites its node numbering: anything already
    /// derived from the element numbering -- cached geometry, per-element tables, a partially
    /// assembled tangent -- describes the old order afterwards and must be rebuilt. Renumber
    /// first, derive second. Registered sidesets and edgesets ARE remapped, and any extra sidesets
    /// passed in `sidesets` are remapped with them.
    ///
    /// PRECONDITION ON THE SWEEP. A colour is conflict-free within its block. Blocks that share
    /// nodes must not be swept concurrently; sweeping them one after another, each with its own
    /// colour loop, is safe and is what the callers here do.
    ///
    /// ORDER MATTERS ON THE WAY IN, TOO. Run this after any mesh reordering, not before: `SFC`
    /// permutes elements and nodes both, and the colouring's quality depends on the element order
    /// it is given (see the note on the visit order in the implementation).
    class ElementColoring final {
    public:
        ElementColoring();
        ~ElementColoring();

        /// `block_names` empty means every block. `sidesets` are remapped along with the
        /// registered ones when `modify_mesh` is set, and ignored otherwise.
        static std::shared_ptr<ElementColoring> create(
                const std::shared_ptr<Mesh>                  &mesh,
                const std::vector<std::string>               &block_names = {},
                const bool                                    modify_mesh = false,
                const std::vector<std::shared_ptr<Sideset>>  &sidesets    = {});

        std::shared_ptr<Mesh> mesh() const;

        ptrdiff_t   n_blocks() const;
        std::string block_name(const int block_idx) const;

        int n_colors(const int block_idx) const;

        /// One colour per element of the block, in the element numbering the mesh has NOW -- so
        /// after a `modify_mesh` build this is sorted, and `color_ptr` is the more useful view.
        SharedBuffer<idx_t> colors(const int block_idx) const;

        /// Colour `c` spans elements `[color_ptr[c], color_ptr[c + 1])`. Only meaningful after a
        /// `modify_mesh` build; without one the elements of a colour are scattered and this
        /// returns null.
        SharedBuffer<ptrdiff_t> color_ptr(const int block_idx) const;

        /// The spread between the smallest and largest colour. Every colour is a
        /// barrier-terminated parallel region, so a colour far below the others is a barrier paid
        /// for a partly idle sweep, and this is what says whether that is happening.
        ptrdiff_t min_elements_per_color(const int block_idx) const;
        ptrdiff_t max_elements_per_color(const int block_idx) const;

        void print(std::ostream &os = std::cout, const int verbosity = 0) const;

    private:
        class Block;
        class Impl;
        std::unique_ptr<Impl> impl_;
    };

}  // namespace smesh

#endif  // SMESH_ELEMENT_COLORING_HPP
