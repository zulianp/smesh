#!/bin/bash
#SBATCH --job-name=mpi_init
#SBATCH --account=c40
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=288
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00
#SBATCH --output=slurm-mpi_init-%j.out
#SBATCH --error=slurm-mpi_init-%j.err
#SBATCH --exclusive
#SBATCH --partition=normal

set -euo pipefail

export MPICH_GPU_SUPPORT_ENABLED=0

WORKDIR="${SLURM_SUBMIT_DIR:-$PWD}"
STEM="${WORKDIR}/mpi_init_bench.$$"
SRC="${STEM}.c"
BIN="${STEM}"

cleanup() { rm -f "$SRC" "$BIN"; }
trap cleanup EXIT

cat >"$SRC" <<'EOF'
#include <mpi.h>
#include <stdio.h>
#include <time.h>

static double walltime(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

int main(int argc, char **argv) {
    double t0 = walltime();
    MPI_Init(&argc, &argv);
    double t1 = walltime();

    MPI_Barrier(MPI_COMM_WORLD);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        printf("MPI_Init time = %.6f s\n", t1 - t0);
    }

    MPI_Finalize();
    return 0;
}
EOF

if command -v mpicc >/dev/null 2>&1; then
    mpicc -O2 -o "$BIN" "$SRC"
else
    cc -O2 -o "$BIN" "$SRC"
fi

echo "#---------------#"
date
echo "SLURM_JOB_ID=${SLURM_JOB_ID:-} SLURM_NTASKS=${SLURM_NTASKS:-}"
echo "#---------------#"

srun "$BIN"

echo "#---------------#"
date
echo "#---------------#"
