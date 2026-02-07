#include <stdio.h>
#include "smesh_path.hpp"

using namespace smesh;

int main(int argc, char** argv) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <block_size> <input_folder> <output_folder>\n", argv[0]);
        return 1;
    }

    // int block_size = std::atoi(argv[1]);
    // Path input_folder(argv[2]);
    // Path output_folder(argv[3]);

    return SMESH_SUCCESS;
}