#ifndef SMESH_ARRAY_COUNT_HPP
#define SMESH_ARRAY_COUNT_HPP

#include "smesh_path.hpp"

#include <cstddef>

namespace smesh {

/// Number of values stored in an array file. Raw files: the file size over the size of the type
/// named by the extension (`raw_value_bytes` for the untyped `.raw`). Compressed files
/// (`<stem>.<dtype>.smz`): the count recorded in their header.
/// Returns -1 (with a message on stderr) if the file cannot be read or its type is unknown.
ptrdiff_t array_count_from_file(const Path &path, size_t raw_value_bytes);

} // namespace smesh

#endif // SMESH_ARRAY_COUNT_HPP
