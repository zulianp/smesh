#include "smesh_output.hpp"
#include "smesh_common.hpp"
#include "smesh_env.hpp"
#include "smesh_glob.hpp"
#include "smesh_smz_io.hpp"
#include "smesh_write.hpp"

#if defined(SMESH_ENABLE_MPI)
#include "smesh_distributed_base.hpp"
#include "smesh_distributed_write.hpp"
#include "smesh_smz_mpi.hpp"
#endif

namespace smesh {

class Output::Impl {
public:
  std::shared_ptr<Mesh> mesh;
  Path path;
  compression::CompressionOptions compression;
  Impl(const std::shared_ptr<Mesh> &mesh, const Path &path)
      : mesh(mesh), path(path),
        compression(compression::CompressionOptions::from_env("SMESH_OUTPUT_COMPRESS")) {}
  ~Impl() {}
};

namespace {
/// Serial compressed write of n values of `type`.
int smz_write_typed(const Path &path, const enum PrimitiveType type, const void *const data,
                    const ptrdiff_t n, const compression::CompressionOptions &opts) {
  using namespace compression;
  switch (type) {
  case SMESH_FLOAT16: return smz_write<f16>(path, static_cast<const f16 *>(data), n, opts);
  case SMESH_FLOAT32: return smz_write<f32>(path, static_cast<const f32 *>(data), n, opts);
  case SMESH_FLOAT64: return smz_write<f64>(path, static_cast<const f64 *>(data), n, opts);
  case SMESH_INT8: return smz_write<i8>(path, static_cast<const i8 *>(data), n, opts);
  case SMESH_INT16: return smz_write<i16>(path, static_cast<const i16 *>(data), n, opts);
  case SMESH_INT32: return smz_write<i32>(path, static_cast<const i32 *>(data), n, opts);
  case SMESH_INT64: return smz_write<i64>(path, static_cast<const i64 *>(data), n, opts);
  case SMESH_UINT8: return smz_write<u8>(path, static_cast<const u8 *>(data), n, opts);
  case SMESH_UINT16: return smz_write<u16>(path, static_cast<const u16 *>(data), n, opts);
  case SMESH_UINT32: return smz_write<u32>(path, static_cast<const u32 *>(data), n, opts);
  case SMESH_UINT64: return smz_write<u64>(path, static_cast<const u64 *>(data), n, opts);
  default: break;
  }
  fprintf(stderr, "Output: type %s cannot be compressed\n", str(type).c_str());
  return SMESH_FAILURE;
}
} // namespace

std::shared_ptr<Mesh> Output::mesh() const { return impl_->mesh; }

Output::Output(const std::shared_ptr<Mesh> &mesh, const Path &path)
    : impl_(std::make_unique<Impl>(mesh, path)) {}

Output::~Output() = default;

std::shared_ptr<Output> Output::create(const std::shared_ptr<Mesh> &mesh,
                                       const Path &path) {
  return std::make_shared<Output>(mesh, path);
}

void Output::set_compression(const compression::CompressionOptions &opts) {
  impl_->compression = opts;
}

const compression::CompressionOptions &Output::compression() const {
  return impl_->compression;
}

int Output::write_nodal(const std::string &field_name,
                        const enum PrimitiveType type, const void *const data,
                        const int block_size) {
  if (block_size != 1) {
    SMESH_ERROR("write_nodal: Block size must be 1 for nodal data");
    return SMESH_FAILURE;
  }
  return write_nodal(field_name, type, data, impl_->compression);
}

int Output::write_elemental(const std::string &field_name,
                            const enum PrimitiveType type,
                            const void *const data, const int block_size) {
  if (block_size != 1) {
    SMESH_ERROR("write_elemental: Block size must be 1 for elemental data");
    return SMESH_FAILURE;
  }
  return write_elemental(field_name, type, data, impl_->compression);
}

int Output::write_nodal(const std::string &field_name,
                        const enum PrimitiveType type, const void *const data,
                        const compression::CompressionOptions &opts) {
  auto mesh = this->mesh();
  if (mesh->comm()->rank() == 0) {
    create_directory(impl_->path);
  }
  mesh->comm()->barrier();

  auto path = Path(impl_->path) /
              (field_name + "." + str(type) + compression::compressed_suffix(&opts));

#if defined(SMESH_ENABLE_MPI)
  if (mesh->comm()->size() > 1) {
    auto dist = mesh->distributed();
    auto data_type = mpi_type_from_primitive_type(type);
    if (opts.enabled) {
      return compression::smz_write_mapped_mpi_type<large_idx_t>(
          mesh->comm()->get(), path, dist->n_nodes_owned(), dist->n_nodes_global(),
          dist->node_mapping()->data(), data_type, data, opts);
    }
    return write_mapped_field(
      mesh->comm()->get(), path, dist->n_nodes_owned(), dist->n_nodes_global(),
      dist->node_mapping()->data(), data_type, data);
  }
#endif
  if (opts.enabled) {
    return smz_write_typed(path, type, data, mesh->n_nodes(), opts);
  }
  return array_write(path, type, data, mesh->n_nodes());
}

int Output::write_elemental(const std::string &field_name,
                            const enum PrimitiveType type,
                            const void *const data,
                            const compression::CompressionOptions &opts) {
  auto mesh = this->mesh();
  if (mesh->comm()->rank() == 0) {
    create_directory(impl_->path);
  }
  mesh->comm()->barrier();

  auto path = Path(impl_->path) /
              (field_name + "." + str(type) + compression::compressed_suffix(&opts));

#if defined(SMESH_ENABLE_MPI)
  if (mesh->comm()->size() > 1) {
    auto dist = mesh->distributed();
    auto data_type = mpi_type_from_primitive_type(type);
    if (opts.enabled) {
      return compression::smz_write_mapped_mpi_type<large_idx_t>(
          mesh->comm()->get(), path, dist->n_elements_owned(), dist->n_elements_global(),
          dist->element_mapping()->data(), data_type, data, opts);
    }
    return write_mapped_field(mesh->comm()->get(), path, dist->n_elements_owned(), dist->n_elements_global(),
                              dist->element_mapping()->data(), data_type, data);
  }
#endif
  if (opts.enabled) {
    return smz_write_typed(path, type, data, mesh->n_elements(), opts);
  }
  return array_write(path, type, data, mesh->n_elements());
}

} // namespace smesh
