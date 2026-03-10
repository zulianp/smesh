#ifndef SMESH_GRID_HPP
#define SMESH_GRID_HPP


#include "smesh_base.hpp"
#include "smesh_communicator.hpp"
#include "smesh_types.hpp"
#include "smesh_buffer.hpp"

#include <cstddef>
#include <memory>
#include <string>

namespace smesh {
    template <class T>
    class Grid;

    template <class T>
    using SharedGrid = std::shared_ptr<Grid<T>>;

    template <class T>
    class Grid {
    public:
        static std::shared_ptr<Grid> create_from_file(const std::shared_ptr<Communicator> &comm, const std::string &path);

        static std::shared_ptr<Grid> create(const std::shared_ptr<Communicator> &comm,
                                            const ptrdiff_t                      nx,
                                            const ptrdiff_t                      ny,
                                            const ptrdiff_t                      nz,
                                            const geom_t                         xmin,
                                            const geom_t                         ymin,
                                            const geom_t                         zmin,
                                            const geom_t                         xmax,
                                            const geom_t                         ymax,
                                            const geom_t                         zmax);

        int to_file(const std::string &folder);

        Grid(const std::shared_ptr<Communicator> &comm);
        ~Grid();

        ptrdiff_t stride(int dim) const;
        ptrdiff_t extent(int dim) const;
        ptrdiff_t size() const;
        int       spatial_dimension() const;
        int       block_size() const;

        const ptrdiff_t *nlocal() const;
        const ptrdiff_t *nglobal() const;
        const ptrdiff_t *stride() const;

        const geom_t *origin() const;
        const geom_t *delta() const;

        SharedBuffer<T> buffer();
        T              *data();

        MPI_Datatype mpi_data_type() const;

#ifdef SMESH_ENABLE_CUDA
        std::shared_ptr<Grid<T>> to_device() const;
#endif

        enum MemorySpace mem_space() const;

        // private:
        class Impl;
        std::unique_ptr<Impl> impl_;
    };

    std::shared_ptr<Grid<geom_t>> create_sdf(const std::shared_ptr<Communicator>                            &comm,
                                             const ptrdiff_t                                                 nx,
                                             const ptrdiff_t                                                 ny,
                                             const ptrdiff_t                                                 nz,
                                             const geom_t                                                    xmin,
                                             const geom_t                                                    ymin,
                                             const geom_t                                                    zmin,
                                             const geom_t                                                    xmax,
                                             const geom_t                                                    ymax,
                                             const geom_t                                                    zmax,
                                             std::function<geom_t(const geom_t, const geom_t, const geom_t)> f);

    template <typename T>
    SharedGrid<T> to_device(const SharedGrid<T> &grid) {
#ifdef SMESH_ENABLE_CUDA
        if (grid->mem_space() == MEMORY_SPACE_DEVICE) {
            return grid->to_device();
        }
#endif
        return grid;
    }

}  // namespace smesh

#endif  // SMESH_GRID_HPP
