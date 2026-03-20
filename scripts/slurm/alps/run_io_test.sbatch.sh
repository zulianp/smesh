#!/bin/bash
#SBATCH --job-name=io_test
#SBATCH --account=c40
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=288
#SBATCH --cpus-per-task=1
#SBATCH --time=00:20:00
#SBATCH --output=slurm-io_test-%j.out
#SBATCH --error=slurm-io_test-%j.err
#SBATCH --exclusive
#SBATCH --partition=normal

set -euo pipefail

export MPICH_GPU_SUPPORT_ENABLED=0
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=true

echo "#---------------#"
date
echo "#---------------#"

# Compile with MPI and OpenMP for hybrid parallelism
srun ./io_test /capstor/scratch/cscs/zulianp/meshes/tet4_cube_768000000 tet4_cube_768000000_data

echo "#---------------#"
date
echo "#---------------#"