#ifndef SMESH_GLOB_HPP
#define SMESH_GLOB_HPP

#include <string>
#include <vector>
#include <string_view>

#include "smesh_path.hpp"

namespace smesh {
	std::vector<std::string> find_files(const std::string_view &pattern);
	size_t count_files(const std::string_view &pattern);
	int create_directory(const std::string_view &path);
	int create_directory(const Path &path);
}

#endif //SMESH_GLOB_HPP
