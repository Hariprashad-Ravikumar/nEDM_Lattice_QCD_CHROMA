#!/bin/bash
#    Begin SBATCH directives
#SBATCH -A m2322_g
#SBATCH -o ./slurm_log.lss
#SBATCH -e ./slurm_err.lss
#SBATCH --nodes=8
#SBATCH -t 00:15:00
#SBATCH --qos=debug
##SBATCH -q regular
#SBATCH --exclusive
#SBATCH -C gpu
#SBATCH --gpus-per-node=4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hari1729@nmsu.edu


source /pscratch/sd/h/hari_8/nEDM_project_LANL/build_script_28Feb2025/env.sh
module load PrgEnv-gnu
###module load cudatoolkit/11.7
module load craype-accel-nvidia80
module load cudatoolkit/12.2

export CRAY_ACCEL_TARGET=nvidia80
export CRAY_CPU_TARGET=x86-64
export MPICH_GPU_SUPPORT_ENABLED=1


run_job(){




ppn=4
let processes=nodes*ppn
let nthread=64/ppn
BIN=/pscratch/sd/h/hari_8/nEDM_project_LANL/build_script_28Feb2025/install/chroma_lanl/bin/chroma
#srun -N 8 -n 32 -c 16 --ntasks-per-node=4 $BIN -geom 1 2 2 8 -i input_prop2ptGFLOW.xml  > logGflowprop2pt3pt
srun -N 8 -n 32 -c 16 --ntasks-per-node=4 $BIN -geom 1 2 2 8 -i input_propBUILDINGBLOCK.xml  > log2pt3pt
} 
run_job
