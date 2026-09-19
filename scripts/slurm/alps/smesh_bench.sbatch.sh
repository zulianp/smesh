#!/bin/bash
#SBATCH --job-name=smesh_bench
#SBATCH --account=c40
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=72
#SBATCH --time=02:00:00
#SBATCH --output=slurm-smesh_bench-%j.out
#SBATCH --error=slurm-smesh_bench-%j.err
#SBATCH --exclusive
#SBATCH --partition=normal

set -euo pipefail

export MPICH_GPU_SUPPORT_ENABLED=0
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-72}
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export SMESH_REORDER=0
export SMESH_CREATE_SIDESETS=0
export SMESH_BENCH_REPEAT=${SMESH_BENCH_REPEAT:-3}
export SMESH_REFINEMENT_LEVELS=${SMESH_REFINEMENT_LEVELS:-1}

BENCH_BIN=${SMESH_BENCH_BIN:-${SMESH_PREFIX:+$SMESH_PREFIX/bin/}smesh_bench}
if [[ ! -x "$BENCH_BIN" ]]; then
  if command -v smesh_bench >/dev/null 2>&1; then
    BENCH_BIN=$(command -v smesh_bench)
  elif [[ -x ./smesh_bench ]]; then
    BENCH_BIN=./smesh_bench
  fi
fi

HEX_N=${SMESH_BENCH_N:-512}
TET_N=${SMESH_BENCH_TET_N:-256}

SCRATCH_ROOT=${SMESH_BENCH_SCRATCH:-${SCRATCH:-/capstor/scratch/cscs/zulianp}/smesh_bench}
WORKDIR="${SCRATCH_ROOT}/${SLURM_JOB_ID:-local}"
mkdir -p "${WORKDIR}"

SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
CSV_DIR="${SMESH_BENCH_CSV_DIR:-${SUBMIT_DIR}/docs/bench}"
mkdir -p "${CSV_DIR}"
CSV="${SMESH_BENCH_CSV:-${CSV_DIR}/results.${SLURM_JOB_ID:-manual}.csv}"

HEX_MESH="${WORKDIR}/hex8_${HEX_N}"
TET_MESH="${WORKDIR}/tet4_${TET_N}"

echo "#---------------#"
date
echo "SLURM_JOB_ID=${SLURM_JOB_ID:-}"
echo "nodes=${SLURM_JOB_NUM_NODES:-} ntasks=${SLURM_NTASKS:-} ntasks_per_node=${SLURM_NTASKS_PER_NODE:-}"
echo "cpus_per_task=${SLURM_CPUS_PER_TASK:-} OMP_NUM_THREADS=${OMP_NUM_THREADS}"
echo "bin=${BENCH_BIN}"
echo "csv=${CSV}"
echo "hex_n=${HEX_N} tet_n=${TET_N}"
echo "workdir=${WORKDIR}"
echo "#---------------#"

run() {
  echo "# $* "
  date
  srun "${BENCH_BIN}" "$@"
}

run generate HEX8 "${HEX_N}" "${CSV}"
run io HEX8 "${HEX_N}" "${CSV}" "${HEX_MESH}"
run refine HEX8 "${HEX_N}" "${CSV}"
run generate TET4 "${TET_N}" "${CSV}"
run io TET4 "${TET_N}" "${CSV}" "${TET_MESH}"
run promote TET4 "${TET_N}" "${CSV}"
run refine TET4 "${TET_N}" "${CSV}"

echo "#---------------#"
date
echo "#---------------#"
