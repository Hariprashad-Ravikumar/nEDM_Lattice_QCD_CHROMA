#!/bin/bash
#    Begin SBATCH directives
#SBATCH -A m2322_g
#SBATCH -o ./slurm_log.lss
#SBATCH -e ./slurm_err.lss
#SBATCH --nodes=3
#SBATCH -t 10:00:00
#SBATCH --qos=regular
##SBATCH -q regular
#SBATCH --mail-user=hari1729@nmsu.edu
#SBATCH --exclusive
#SBATCH -C gpu
#SBATCH --gpus-per-node=4

source /pscratch/sd/h/hari_8/nEDM_project_LANL/build_script_28Feb2025/env.sh
module load PrgEnv-gnu
###module load cudatoolkit/11.7
module load craype-accel-nvidia80
module load cudatoolkit/12.2

export CRAY_ACCEL_TARGET=nvidia80
export CRAY_CPU_TARGET=x86-64
export MPICH_GPU_SUPPORT_ENABLED=1


run_job(){



geom="3 1 2 2"
ppn=4
let processes=nodes*ppn
let nthread=64/ppn
#BIN=/global/homes/f/fche123/fche123/package-20231108_cpu/cpu-Milan_gcc/install/chroma_lanl-double/bin/chroma
BIN=/pscratch/sd/h/hari_8/nEDM_project_LANL/build_script_28Feb2025/install/chroma_lanl/bin/chroma
#BIN=/global/common/software/mp7/USQCD/scripts_Perlmutter_CPU_GNU_newflags/build/chroma/bin/chroma
#srun -n $processes -ppn $ppn $BIN --mpi ${geomm} -nthread $nthread -i input.a_2004.xml -o out.a_2004.xml > log111
#srun -N 2 -n 32 --cpu_bind=core $BIN  -i input.a_2004.xml -o out.a_2004.xml -x 64 -y 64 -z 64 -t 96 -by 8 -bz 8 -pxy 1 -pxyz 0 -c 2 -sy 1 -sz 1 -minct 2 > log111
srun -N 3 -n 12 $BIN -geom 3 1 2 2 -i input_propagator.xml > logforprop2pt
#srun -N 4  -n 512 -c 2 --cpu_bind=cores $BIN -geom 4 4 4 8 -i input.a_2004.xml > log111
} 
run_job
