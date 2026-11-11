#ifndef SMESH_OUTPUT_HPP
#define SMESH_OUTPUT_HPP

#include "smesh_base.hpp"
#include "smesh_buffer.hpp"
#include "smesh_forward_declarations.hpp"
#include "smesh_types.hpp"
#include "smesh_mesh.hpp"
#include "smesh_compression_options.hpp"

#include <memory>

namespace smesh {
class Output {
public:
  Output(const std::shared_ptr<Mesh> &mesh, const Path &path);
  ~Output();

  static std::shared_ptr<Output> create(const std::shared_ptr<Mesh> &mesh,
                                        const Path &path);

  std::shared_ptr<Mesh> mesh() const;

  /// Compression of the fields written from now on (default: SMESH_OUTPUT_COMPRESS, off). Files are
  /// then `<name>.<dtype>.smz`; e.g. "range:1e-4" bounds the error by 1e-4 of each field's range.
  void                                   set_compression(const compression::CompressionOptions &opts);
  const compression::CompressionOptions &compression() const;

  int write_nodal(const std::string &field_name, const enum PrimitiveType type, const void *const data,
                  const int block_size = 1);

  int write_elemental(const std::string &field_name, const enum PrimitiveType type, const void *const data,
                      const int block_size = 1);

  /// As above, with explicit compression for this field.
  int write_nodal(const std::string &field_name, const enum PrimitiveType type, const void *const data,
                  const compression::CompressionOptions &opts);
  int write_elemental(const std::string &field_name, const enum PrimitiveType type, const void *const data,
                      const compression::CompressionOptions &opts);

  template <typename T> int write_nodal(const std::string &field_name, const SharedBuffer<T> &data) {
    ptrdiff_t n = data->size();
    int block_size = n / mesh()->n_nodes();
    return write_nodal(field_name, TypeToEnum<T>::value(), data->data(), block_size);
  }

  template <typename T> int write_elemental(const std::string &field_name, const SharedBuffer<T> &data) {
    ptrdiff_t n = data->size();
    int block_size = n / mesh()->n_elements();
    return write_elemental(field_name, TypeToEnum<T>::value(), data->data(), block_size);
  }

  class Impl;

private:
  std::unique_ptr<Impl> impl_;
};

} // namespace smesh

#endif // SMESH_OUTPUT_HPP
