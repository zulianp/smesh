#!/bin/bash
#SBATCH --job-name=mpi_io_bench
#SBATCH --account=c40
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=288
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --output=slurm-mpi_io_bench-%j.out
#SBATCH --error=slurm-mpi_io_bench-%j.err
#SBATCH --exclusive
#SBATCH --partition=normal

set -euo pipefail

export MPICH_GPU_SUPPORT_ENABLED=0
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OMP_PROC_BIND=true

BENCH_BIN=${BENCH_BIN:-./mpi_io_bench}
INPUT=${INPUT:?set INPUT=/path/to/array.raw}
DTYPE=${DTYPE:-float32}
SEGMENT_MB_LIST=${SEGMENT_MB_LIST:-256 512 1024}
LEADERS_PER_NODE_LIST=${LEADERS_PER_NODE_LIST:-1 2 4}

dtype_size() {
  case "$1" in
    float32|f32|int32|i32|uint32|u32) echo 4 ;;
    float64|f64|int64|i64|uint64|u64) echo 8 ;;
    uint16|u16) echo 2 ;;
    uint8|u8) echo 1 ;;
    *)
      echo "Unsupported DTYPE: $1" >&2
      exit 1
      ;;
  esac
}

TYPE_SIZE=$(dtype_size "${DTYPE}")

segment_elems() {
  local mb=$1
  echo $(( mb * 1024 * 1024 / TYPE_SIZE ))
}

run_case() {
  local strategy=$1
  local segment=${2:-2147483647}
  local leaders=${3:-1}

  echo "# strategy=${strategy} segment_elems=${segment} leaders_per_node=${leaders}"
  date
  srun "${BENCH_BIN}" "${INPUT}" "${DTYPE}" "${strategy}" "${segment}" "${leaders}"
}

echo "#---------------#"
echo "# mpi_io_bench"
echo "# bin=${BENCH_BIN}"
echo "# input=${INPUT}"
echo "# dtype=${DTYPE}"
echo "# nodes=${SLURM_JOB_NUM_NODES:-unknown}"
echo "# ntasks=${SLURM_NTASKS:-unknown}"
echo "# ntasks_per_node=${SLURM_NTASKS_PER_NODE:-unknown}"
date
echo "#---------------#"

run_case matrixio

for segment_mb in ${SEGMENT_MB_LIST}; do
  segment=$(segment_elems "${segment_mb}")
  run_case matrixio_segmented "${segment}"
  run_case mpi_collective "${segment}"
  run_case mpi_independent "${segment}"

  for leaders in ${LEADERS_PER_NODE_LIST}; do
    run_case node_aggregated "${segment}" "${leaders}"
  done
done

echo "#---------------#"
date
echo "#---------------#"
