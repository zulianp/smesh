#!/bin/bash
#SBATCH --job-name=gencube
#SBATCH --account=c40
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=72
#SBATCH --time=04:00:00
#SBATCH --output=slurm-gencube-%j.out
#SBATCH --error=slurm-gencube-%j.err
#SBATCH --exclusive
#SBATCH --partition=normal

set -euo pipefail

export MPICH_GPU_SUPPORT_ENABLED=1
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=true

export PATH=/capstor/scratch/cscs/zulianp/installations/smesh/bin:$PATH

echo "#---------------#"
date
echo "#---------------#"

export SMESH_REORDER=0
export SMESH_CREATE_SIDESETS=0

set -x

./cube TET4 400 400 400    -10 -10 -10 10 10 10 tet4_cube_768000000
./cube TET4 1000 1000 1000 -10 -10 -10 10 10 10 tet4_cube_12000000000
./cube HEX8 5000 5000 5000 -10 -10 -10 10 10 10 hex8_cube_125000000000

set +x

echo "#---------------#"
date
echo "#---------------#"