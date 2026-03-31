# Test case 1

To create test 1 generate a mesh with 12 Giga elements with the following command (you need 381GB of RAM as well as disk space)


```bash
 ./cube TET4 1000 1000 1000 -10 -10 -10 10 10 10 tet4_cube
```

```yaml
# smesh mesh meta file
spatial_dimension: 3
elem_num_nodes: 4
element_type: TET4
n_elements: 12000000000
n_nodes: 2003003001
elements:
- i0: i0.int64
- i1: i1.int64
- i2: i2.int64
- i3: i3.int64
points:
- x: x.float32
- y: y.float32
- z: z.float32
rpath: true
```

One node per socket slurm script

```bash
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
```

## Breakdown

The executable takes `296.793 [s]` on the Alps supercomputer (8 nodes, 32 MPI ranks).
MPI initialization/finalization takes `5.57 [s]` which is `1.8` percent of runtime.

MPI API (reads) 55.6 percent of time:
- Read x,y,z files `16.0297 [s]` 
- Read indices files `149.149 [s]`

all_to_allv_64 routines are accounting for approximately `20.2 [s]`

MPI-based routines account approximately for 64 percent of set-up time
`(5.57 + 16.0297 + 149.149 + 20.2) [s] = 190.9 [s]`





## Options


