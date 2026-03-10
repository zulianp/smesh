
#include "smesh_mask.hpp"
#include "smesh_context.hpp"
#include "smesh_types.hpp"
#include "smesh_buffer.hpp"
#include "smesh_write.hpp"
#include "smesh_env.hpp"

#include <stdio.h>

using namespace smesh;

int main(int argc, char *argv[]) {
    auto ctx = initialize_serial(argc, argv);

    if (argc != 4) {
        fprintf(stderr, "usage: %s <n> <input> <output>\n", argv[0]);
        return EXIT_FAILURE;
    }

    bool SFEM_SKIP_WRITE = Env::read("SFEM_SKIP_WRITE", false);
    const ptrdiff_t n = atol(argv[1]);
    auto path_input = Path(argv[2]);
    auto path_output = Path(argv[3]);

    auto idx = Buffer<idx_t>::from_file(path_input);
    auto data = idx->data();
    const ptrdiff_t size = idx->size();

    auto mask = create_host_buffer<mask_t>(mask_count(n));
    auto m = mask->data();

    double tick = time_seconds();

#pragma omp parallel for
    for (ptrdiff_t i = 0; i < size; i++) {
        mask_set(data[i], m);
    }

    double tock = time_seconds();
    double elapsed = tock - tick;

    ptrdiff_t output_bytes = mask_count(n) * sizeof(mask_t);
    ptrdiff_t intput_bytes = size * sizeof(idx_t);
    printf("Mask set TTS %g [s] BW %g [GB/s]\n",
           elapsed,
           1e-9 * (intput_bytes + output_bytes) / elapsed);

#ifndef NDEBUG
    for (ptrdiff_t i = 0; i < size; i++) {
        int must_be_true = mask_get(data[i], m);
        SMESH_ASSERT(must_be_true);
    }
#endif

    if (!SFEM_SKIP_WRITE) {
        mask->to_file(path_output);
    }

#ifndef NDEBUG
    for (ptrdiff_t i = 0; i < size; i++) {
        mask_unset(data[i], m);
        int must_be_false = !mask_get(data[i], m);
        SMESH_ASSERT(must_be_false);
    }
#endif

    return SMESH_SUCCESS;
}
