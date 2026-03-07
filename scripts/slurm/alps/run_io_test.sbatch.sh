#!/bin/bash
#SBATCH --job-name=io_test
#SBATCH --account=c40
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=72
#SBATCH --time=00:20:00
#SBATCH --output=slurm-io_test-%j.out
#SBATCH --error=slurm-io_test-%j.err
#SBATCH --exclusive
#SBATCH --partition=normal

set -euo pipefail

export MPICH_GPU_SUPPORT_ENABLED=1
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=true

echo "#---------------#"
date
echo "#---------------#"

# Compile with MPI and OpenMP for hybrid parallelism
srun ./io_test tet4_cube tet4_cube_data

echo "#---------------#"
date
echo "#---------------#"