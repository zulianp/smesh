#include "smesh_grid.hpp"
#include "smesh_device_buffer.hpp"
#include "smesh_glob.hpp"
#include "smesh_tracer.hpp"

#ifdef SMESH_ENABLE_MPI
#include "matrixio_ndarray.h"
#include "smesh_distributed_base.hpp"
#endif

#ifdef SMESH_ENABLE_RYAML
#include <ryml.hpp>
#include <ryml_std.hpp>
#endif

#include <fstream>
#include <iterator>
#include <sstream>

namespace smesh {

template <class T> class Grid<T>::Impl {
public:
  std::shared_ptr<Communicator> comm;
  std::shared_ptr<Buffer<T>> field;
  int block_size{1};
  int spatial_dimension{3};

  ptrdiff_t nlocal[3];
  ptrdiff_t nglobal[3];
  ptrdiff_t stride[3];

  // Grid geometry
  geom_t origin[3];
  geom_t delta[3];

  enum MemorySpace mem_space { MEMORY_SPACE_HOST };
};

#ifdef SMESH_ENABLE_MPI
template <class T> MPI_Datatype Grid<T>::mpi_data_type() const {
  return mpi_type<T>();
}
#endif

template <class T>
Grid<T>::Grid(const std::shared_ptr<Communicator> &comm)
    : impl_(std::make_unique<Impl>()) {
  impl_->comm = comm;
}

template <class T> Grid<T>::~Grid() = default;

template <class T>
std::shared_ptr<Grid<T>>
Grid<T>::create_from_file(const std::shared_ptr<Communicator> &comm,
                          const std::string &path) {
  auto ret = std::make_unique<Grid<T>>(comm);
  ret->impl_->nlocal[0] = ret->impl_->nlocal[1] = ret->impl_->nlocal[2] = 1;
  ret->impl_->nglobal[0] = ret->impl_->nglobal[1] = ret->impl_->nglobal[2] = 1;
  ret->impl_->origin[0] = ret->impl_->origin[1] = ret->impl_->origin[2] = 0;
  ret->impl_->delta[0] = ret->impl_->delta[1] = ret->impl_->delta[2] = 1;

  std::string field_path;
#if defined(SMESH_ENABLE_RYAML)
  {
    std::ifstream ifs(path + "/meta.yaml", std::ios::binary);
    if (!ifs.good()) {
      SMESH_ERROR("Grid::create_from_file: Unable to read %s\n",
                  (path + "/meta.yaml").c_str());
      return nullptr;
    }

    std::string yaml((std::istreambuf_iterator<char>(ifs)),
                     std::istreambuf_iterator<char>());
    if (yaml.empty()) {
      SMESH_ERROR("Grid::create_from_file: Empty metadata file %s\n",
                  (path + "/meta.yaml").c_str());
      return nullptr;
    }

    auto tree = ryml::parse_in_arena(ryml::to_csubstr(yaml));
    auto root = tree.rootref();
    auto require = [&root](const char *key, auto &value) {
      if (!root.has_child(key)) {
        return false;
      }

      root[key] >> value;
      return true;
    };

    auto get = [&root](const char *key, auto &value) {
      if (!root.has_child(key)) {
        return false;
      }

      root[key] >> value;
      return true;
    };

    if (!require("path", field_path) ||
        !require("spatial_dimension", ret->impl_->spatial_dimension) ||
        !require("nx", ret->impl_->nglobal[0]) ||
        !require("ox", ret->impl_->origin[0]) ||
        !require("dx", ret->impl_->delta[0])) {
      SMESH_ERROR("Grid::create_from_file: Missing required metadata in %s\n",
                  (path + "/meta.yaml").c_str());
      return nullptr;
    }

    if (ret->impl_->spatial_dimension > 1 &&
        (!require("ny", ret->impl_->nglobal[1]) ||
         !require("oy", ret->impl_->origin[1]) ||
         !require("dy", ret->impl_->delta[1]))) {
      SMESH_ERROR("Grid::create_from_file: Missing required metadata in %s\n",
                  (path + "/meta.yaml").c_str());
      return nullptr;
    }

    if (ret->impl_->spatial_dimension > 2 &&
        (!require("nz", ret->impl_->nglobal[2]) ||
         !require("oz", ret->impl_->origin[2]) ||
         !require("dz", ret->impl_->delta[2]))) {
      SMESH_ERROR("Grid::create_from_file: Missing required metadata in %s\n",
                  (path + "/meta.yaml").c_str());
      return nullptr;
    }

    get("block_size", ret->impl_->block_size);
  }
#else
  {
    std::ifstream ifs(path + "/meta.yaml");
    if (!ifs.good()) {
      SMESH_ERROR("Grid::create_from_file: Unable to read %s\n",
                  (path + "/meta.yaml").c_str());
      return nullptr;
    }

    auto read_scalar = [](const std::string &line, const char *key,
                          auto &value) -> bool {
      if (line.rfind(key, 0) != 0) {
        return false;
      }

      std::istringstream is(line.substr(std::char_traits<char>::length(key)));
      is >> value;
      return !is.fail();
    };

    std::string line;
    while (std::getline(ifs, line)) {
      if (line.rfind("path:", 0) == 0) {
        field_path = line.substr(5);
        const auto first = field_path.find_first_not_of(" \t");
        field_path = first == std::string::npos ? std::string()
                                                : field_path.substr(first);
        continue;
      }

      read_scalar(line, "spatial_dimension:", ret->impl_->spatial_dimension) ||
          read_scalar(line, "nx:", ret->impl_->nglobal[0]) ||
          read_scalar(line, "ox:", ret->impl_->origin[0]) ||
          read_scalar(line, "dx:", ret->impl_->delta[0]) ||
          read_scalar(line, "ny:", ret->impl_->nglobal[1]) ||
          read_scalar(line, "oy:", ret->impl_->origin[1]) ||
          read_scalar(line, "dy:", ret->impl_->delta[1]) ||
          read_scalar(line, "nz:", ret->impl_->nglobal[2]) ||
          read_scalar(line, "oz:", ret->impl_->origin[2]) ||
          read_scalar(line, "dz:", ret->impl_->delta[2]) ||
          read_scalar(line, "block_size:", ret->impl_->block_size);
    }

    if (field_path.empty() || ret->impl_->spatial_dimension <= 0) {
      SMESH_ERROR("Grid::create_from_file: Missing required metadata in %s\n",
                  (path + "/meta.yaml").c_str());
      return nullptr;
    }
  }
#endif

  assert(ret->impl_->spatial_dimension > 0);
  assert(ret->impl_->spatial_dimension <= 3);

  ret->impl_->nlocal[0] = ret->impl_->nglobal[0];
  ret->impl_->nlocal[1] = ret->impl_->nglobal[1];
  ret->impl_->nlocal[2] = ret->impl_->nglobal[2];

  if (!field_path.empty() && field_path[0] != '/') {
    field_path = (smesh::Path(path) / field_path).to_string();
  }

#ifdef SMESH_ENABLE_MPI
  {
    T *data = nullptr;
    if (ndarray_create_from_file(
            ret->impl_->comm->get(), field_path.c_str(), ret->mpi_data_type(),
            ret->impl_->spatial_dimension, (void **)&data, ret->impl_->nlocal,
            ret->impl_->nglobal) != SMESH_SUCCESS) {
      SMESH_ERROR("Grid::create_from_file: Unable to read %s\n",
                  field_path.c_str());
      return nullptr;
    }

    ptrdiff_t size = ret->impl_->nlocal[0];
    if (ret->impl_->spatial_dimension > 1)
      size *= ret->impl_->nlocal[1];
    if (ret->impl_->spatial_dimension > 2)
      size *= ret->impl_->nlocal[2];
    ret->impl_->field = Buffer<T>::own(size, data, free, MEMORY_SPACE_HOST);
  }
#else
  ret->impl_->field = Buffer<T>::from_file(smesh::Path(field_path));
#endif

  ret->impl_->stride[0] = 1;
  ret->impl_->stride[1] = ret->impl_->nlocal[0];
  ret->impl_->stride[2] = ret->impl_->nlocal[0] * ret->impl_->nlocal[1];
  ptrdiff_t size = ret->impl_->stride[2] * ret->impl_->nlocal[2];
  T *data = ret->impl_->field->data();

  geom_t SMESH_GRID_SHIFT = 0;
  geom_t SMESH_GRID_SCALE = 1;

  SMESH_READ_ENV(SMESH_GRID_SHIFT, atof);
  SMESH_READ_ENV(SMESH_GRID_SCALE, atof);

#pragma omp parallel for
  for (ptrdiff_t i = 0; i < size; i++) {
    data[i] = SMESH_GRID_SCALE * (data[i] + SMESH_GRID_SHIFT);
  }
  return ret;
}

template <class T> int Grid<T>::to_file(const std::string &folder) {
  std::stringstream ss;

  std::string field_path;
  const std::string field_ext = std::string(smesh::TypeToString<T>::value());
  ss << "path: "
     << "sdf." << field_ext << "\n";
  ss << "spatial_dimension: " << impl_->spatial_dimension << "\n";
  ss << "nx: " << impl_->nglobal[0] << "\n";
  ss << "ox: " << impl_->origin[0] << "\n";
  ss << "dx: " << impl_->delta[0] << "\n";

  if (impl_->spatial_dimension > 1) {
    ss << "ny: " << impl_->nglobal[1] << "\n";
    ss << "oy: " << impl_->origin[1] << "\n";
    ss << "dy: " << impl_->delta[1] << "\n";
  }

  if (impl_->spatial_dimension > 2) {
    ss << "nz: " << impl_->nglobal[2] << "\n";
    ss << "oz: " << impl_->origin[2] << "\n";
    ss << "dz: " << impl_->delta[2] << "\n";
  }

  ss << "block_size: " << impl_->block_size << "\n";

  smesh::create_directory(folder.c_str());
  std::ofstream os(folder + "/meta.yaml");
  if (!os.good())
    return SMESH_FAILURE;
  os << ss.str();
  os.close();

  return impl_->field->to_file(smesh::Path(folder + "/sdf." + field_ext));
}

template <class T> ptrdiff_t Grid<T>::stride(int dim) const {
  return impl_->stride[dim];
}

template <class T> ptrdiff_t Grid<T>::extent(int dim) const {
  return impl_->nlocal[dim];
}

template <class T> ptrdiff_t Grid<T>::size() const {
  return impl_->field->size();
}

template <class T> int Grid<T>::block_size() const { return impl_->block_size; }

template <class T> const ptrdiff_t *Grid<T>::nlocal() const {
  return impl_->nlocal;
}

template <class T> const ptrdiff_t *Grid<T>::nglobal() const {
  return impl_->nglobal;
}

template <class T> const ptrdiff_t *Grid<T>::stride() const {
  return impl_->stride;
}

template <class T> int Grid<T>::spatial_dimension() const {
  return impl_->spatial_dimension;
}

template <class T> std::shared_ptr<Buffer<T>> Grid<T>::buffer() {
  return impl_->field;
}

template <class T> T *Grid<T>::data() { return impl_->field->data(); }

template <class T> const geom_t *Grid<T>::origin() const {
  return impl_->origin;
}

template <class T> const geom_t *Grid<T>::delta() const { return impl_->delta; }

template <class T>
std::shared_ptr<Grid<T>>
Grid<T>::create(const std::shared_ptr<Communicator> &comm, const ptrdiff_t nx,
                const ptrdiff_t ny, const ptrdiff_t nz, const geom_t xmin,
                const geom_t ymin, const geom_t zmin, const geom_t xmax,
                const geom_t ymax, const geom_t zmax) {
  auto ret = std::make_unique<Grid<T>>(comm);

  auto &impl = ret->impl_;

  impl->block_size = 1;
  impl->spatial_dimension = 3;

  // FIXME
  impl->nlocal[0] = nx;
  impl->nlocal[1] = ny;
  impl->nlocal[2] = nz;

  // FIXME
  impl->stride[0] = 1;
  impl->stride[1] = nx;
  impl->stride[2] = nx * ny;

  // FIXME
  impl->nglobal[0] = nx;
  impl->nglobal[1] = ny;
  impl->nglobal[2] = nz;

  // Grid geometry
  impl->origin[0] = xmin;
  impl->origin[1] = ymin;
  impl->origin[2] = zmin;

  impl->delta[0] = (xmax - xmin) / (nx - 1);
  impl->delta[1] = (ymax - ymin) / (ny - 1);
  impl->delta[2] = (zmax - zmin) / (nz - 1);

  impl->field = create_host_buffer<T>(nx * ny * nz);
  return ret;
}

template <typename T> enum MemorySpace Grid<T>::mem_space() const {
  return impl_->mem_space;
}

#ifdef SMESH_ENABLE_CUDA
template <typename T> std::shared_ptr<Grid<T>> Grid<T>::to_device() const {
  auto ret = std::make_shared<Grid<T>>(impl_->comm);

  ret->impl_->field = smesh::to_device(impl_->field);
  ret->impl_->block_size = impl_->block_size;
  ret->impl_->spatial_dimension = impl_->spatial_dimension;

  for (int d = 0; d < 3; d++) {
    ret->impl_->nlocal[d] = impl_->nlocal[d];
    ret->impl_->nglobal[d] = impl_->nglobal[d];
    ret->impl_->stride[d] = impl_->stride[d];
    ret->impl_->origin[d] = impl_->origin[d];
    ret->impl_->delta[d] = impl_->delta[d];
  }

  ret->impl_->mem_space = MEMORY_SPACE_DEVICE;
  return ret;
}
#endif // SMESH_ENABLE_CUDA

std::shared_ptr<Grid<geom_t>>
create_sdf(const std::shared_ptr<Communicator> &comm, const ptrdiff_t nx,
           const ptrdiff_t ny, const ptrdiff_t nz, const geom_t xmin,
           const geom_t ymin, const geom_t zmin, const geom_t xmax,
           const geom_t ymax, const geom_t zmax,
           std::function<geom_t(const geom_t, const geom_t, const geom_t)> f) {
  SMESH_TRACE_SCOPE("Grid::create_sdf");

  auto g = Grid<geom_t>::create(comm, nx, ny, nz, xmin, ymin, zmin, xmax, ymax,
                                zmax);

  const geom_t hx = (xmax - xmin) / (nx - 1);
  const geom_t hy = (ymax - ymin) / (ny - 1);
  const geom_t hz = (zmax - zmin) / (nz - 1);

  geom_t *const field = g->buffer()->data();
  auto stride = g->stride();

#pragma omp parallel for
  for (ptrdiff_t zi = 0; zi < nz; zi++) {
    for (ptrdiff_t yi = 0; yi < ny; yi++) {
      for (ptrdiff_t xi = 0; xi < nx; xi++) {
        const geom_t x = xmin + xi * hx;
        const geom_t y = ymin + yi * hy;
        const geom_t z = zmin + zi * hz;
        const geom_t fx = f(x, y, z);
        field[xi * stride[0] + yi * stride[1] + zi * stride[2]] = fx;
      }
    }
  }

  return g;
}

// Explicit instantiation
template class Grid<float>;
template class Grid<double>;

} // namespace smesh
