#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "smesh_base.hpp"
#include "smesh_context.hpp"
#include "smesh_types.hpp"
#include "smesh_buffer.hpp"
#include "smesh_write.hpp"
#include "smesh_env.hpp"
// #include "smesh_glob.hpp"

using namespace smesh;  

int main(int argc, char *argv[]) {
    auto ctx = initialize_serial(argc, argv);

    if (argc != 4) {
        fprintf(stderr, "usage: %s <alpha> <x> <y>\nApplies y += alpha * x\n", argv[0]);
        return EXIT_FAILURE;
    }

    const real_t alpha = atof(argv[1]);
    auto x_path = Path(argv[2]);
    auto y_path = Path(argv[3]);

    auto x = Buffer<real_t>::from_file(x_path);
    auto y = Buffer<real_t>::from_file(y_path);
    const ptrdiff_t size = y->size();

    auto b_x = x->data();
    auto b_y = y->data();
    for(ptrdiff_t i = 0; i < size; ++i) {
        b_y[i] += alpha * b_x[i];
    }

    y->to_file(y_path);

    return EXIT_SUCCESS;
}
