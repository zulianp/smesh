#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <cassert>
#include <stdexcept>


template <typename T>
void read_array_from_file(const char* fname, std::vector<T>& array) {
   FILE* file = fopen(fname, "rb");
   if (file == NULL) {
    throw std::runtime_error("Failed to open file: " + std::string(fname));
   }

   size_t nbytes = fseek(file, 0, SEEK_END);
   fseek(file, 0, SEEK_SET);
   array.resize(nbytes / sizeof(T));

   assert(nbytes == array.size() * sizeof(T));
   fread(array.data(), sizeof(T), array.size(), file);
   fclose(file);
}

template <typename T>
void write_array_to_file(const char* fname, const std::vector<T>& array) {
    FILE* file = fopen(fname, "wb");
    if (file == NULL) {
        throw std::runtime_error("Failed to open file: " + std::string(fname));
    }
    fwrite(array.data(), sizeof(T), array.size(), file);
    fclose(file);
}

int main(int argc, char** argv) {
    if (argc != 7) {
        std::cerr << "Usage: " << argv[0] << " <i0.int32> <i1.int32> <i2.int32> <x.float32> <y.float32> <z.float32>" << std::endl;
        return 1;
    }

    std::vector<int32_t> i0, i1, i2;
    read_array_from_file(argv[1], i0);
    read_array_from_file(argv[2], i1);
    read_array_from_file(argv[3], i2);

    std::vector<float> x, y, z;
    read_array_from_file(argv[4], x);
    read_array_from_file(argv[5], y);
    read_array_from_file(argv[6], z);

    return 0;
}