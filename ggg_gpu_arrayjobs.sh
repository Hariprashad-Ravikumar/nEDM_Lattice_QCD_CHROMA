#!/bin/bash
#SBATCH -A m2322_g
#SBATCH -o ./slurm_log__gg_%A_%a.lss  # Output file for each array task
#SBATCH -e ./slurm_err__gg_%A_%a.lss  # Error file for each array task
#SBATCH --nodes=8
#SBATCH -t 02:00:00
##SBATCH --qos=debug
#SBATCH -q regular
#SBATCH --exclusive
#SBATCH -C gpu
#SBATCH --gpus-per-node=4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hari1729@nmsu.edu
#SBATCH --array=0-19  # Set the range of tasks based on the number of serial numbers in the text file

# Load necessary modules
source /pscratch/sd/h/hari_8/nEDM_project_LANL/build_script_28Feb2025/env.sh
module load PrgEnv-gnu
module load craype-accel-nvidia80
module load cudatoolkit/12.2

export CRAY_ACCEL_TARGET=nvidia80
export CRAY_CPU_TARGET=x86-64
export MPICH_GPU_SUPPORT_ENABLED=1

# Define the template and file paths
#TEMPLATE=input_cfg_No_HP_GFlow.xml
TEMPLATE=input_ggg_cfg.xml
BIN=/pscratch/sd/h/hari_8/nEDM_project_LANL/build_script_28Feb2025/install/chroma_lanl/bin/chroma

# Path to the file containing the serial numbers
SERIAL_NUMBERS_FILE=serial_numbers.txt

# Read the serial number for this array task
CFG=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" $SERIAL_NUMBERS_FILE)
echo "Running for serial number: ${CFG}"
export CFG
# Create a FIFO file for input substitution
PIPE=./cfg_${CFG}.fifo
rm -f "$PIPE"
mkfifo "$PIPE"
envsubst '${CFG}' < "$TEMPLATE" > "$PIPE" &

# Run the Chroma job
#srun -N 8 -n 32 -c 16 --ntasks-per-node=4 $BIN -geom 1 2 2 8 -i "$PIPE" > log2pt3pt_${CFG}
srun -N 8 -n 32 -c 16 --ntasks-per-node=4 $BIN -geom 1 2 2 8 -i "$PIPE" -o /pscratch/sd/h/hari_8/nEDM_project_LANL/output_3pt2pt_flow/out_qtop_ggg_allLFLOW_cfg_${CFG}.xml > log2pt3pt_${CFG}
# Clean up FIFO after job completes
wait
rm -f "$PIPE"

