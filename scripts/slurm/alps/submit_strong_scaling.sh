#!/usr/bin/env bash
# Submit hybrid strong-scaling jobs: 1 socket (1x72) through 16 nodes (64x72).
# Override sbatch flags on the command line; this script only sets nodes / ntasks-per-node.

set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname "$0")" >/dev/null 2>&1 && pwd)"
SBATCH_SCRIPT="${SCRIPT_DIR}/smesh_bench.sbatch.sh"

ACCOUNT=${ACCOUNT:-c40}
PARTITION=${PARTITION:-normal}
EXTRA_SBATCH_ARGS=${EXTRA_SBATCH_ARGS:-}

submit() {
  local label=$1
  local nodes=$2
  local ntasks_per_node=$3
  echo "submit ${label}: nodes=${nodes} ntasks-per-node=${ntasks_per_node} (cpus-per-task=72)"
  sbatch --job-name="smesh_bench_${label}" \
         --account="${ACCOUNT}" \
         --partition="${PARTITION}" \
         --nodes="${nodes}" \
         --ntasks-per-node="${ntasks_per_node}" \
         --cpus-per-task=72 \
         --exclusive \
         ${EXTRA_SBATCH_ARGS} \
         "${SBATCH_SCRIPT}"
}

# 1 GH200 socket, then 1..16 nodes at 4 ranks/node (one rank per socket).
submit socket 1 1
submit n1 1 4
submit n2 2 4
submit n4 4 4
submit n8 8 4
submit n16 16 4
