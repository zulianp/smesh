#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "smesh_base.hpp"
#include "smesh_context.hpp"
#include "smesh_types.hpp"
#include "smesh_buffer.hpp"
#include "smesh_common.hpp"
#include "smesh_write.hpp"
#include "smesh_glob.hpp"

using namespace smesh;

int main(int argc, char *argv[]) {
    auto ctx = initialize_serial(argc, argv);

    if (argc < 3) {
        fprintf(stderr, "usage: %s <idx.raw> <n> [out=indicator.raw]", argv[0]);
        return EXIT_FAILURE;
    }

    typedef float indicator_t;

    auto output_path = Path("./indicator." + str(TypeToString<indicator_t>::value()));

    if (argc > 3) {
        output_path = Path(argv[3]);
    }

    double tick = time_seconds();

    ptrdiff_t nnodes = atoll(argv[2]);
    auto indicator = create_host_buffer<indicator_t>(nnodes);

    {
        auto b_indicator = indicator->data();
        auto indices = Buffer<idx_t>::from_file(Path(argv[1]));
        auto b_indices = indices->data();

        ptrdiff_t _nope_, ndirichlet;
        for (ptrdiff_t node = 0; node < ndirichlet; ++node) {
            b_indicator[b_indices[node]] = (indicator_t)1;
        }
    }

    indicator->to_file(output_path);

    double tock = time_seconds();
    printf("TTS: %g seconds\n", tock - tick);
    return EXIT_SUCCESS;
}
