#!/bin/bash
#    Begin SBATCH directives
#SBATCH -A m2322_g
#SBATCH -o ./slurm_log4444_3662.lss
#SBATCH -e ./slurm_err4444_3662.lss
#SBATCH --nodes=2
###SBATCH --nodes=1
#SBATCH -t 00:30:00
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
srun -N 2 -n 8 -c 16 --ntasks-per-node=4 $BIN -geom 1 1 2 4 -i input_prop_try3.xml  > logprop_blk4444_3662light
#srun -N 6 -n 24 -c 16 --ntasks-per-node=4 /global/cfs/cdirs/m2322/junsik/qcd_progs/package-20220930_gpu/build-20250124/install/chroma_lanl/bin/chroma -i /global/cfs/cdirs/m2322/junsik/Measurements/NME/a09m130_l6496f211b630m001326m03636m4313/meas-hyp6/input.xml/input.a_14514.LP.05.xml -geom 1 2 2 6 -iogeom 1 1 2 6 >& /global/cfs/cdirs/m2322/junsik/Measurements/NME/a09m130_l6496f211b630m001326m03636m4313/meas-hyp6/job.log/job.log.a_14514.LP.05 &
} 
run_job
