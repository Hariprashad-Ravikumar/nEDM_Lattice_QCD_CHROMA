#!/bin/bash
#    Begin SBATCH directives
#SBATCH -A m2322_g
#SBATCH -o ./slurm_log.lss
#SBATCH -e ./slurm_err.lss
#SBATCH --nodes=8
#SBATCH -t 05:30:00
##SBATCH --qos=debug
#SBATCH -q regular
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
    CFG=$1
    export CFG

    # path to your one‐and‐only template
    TEMPLATE=input_cfg_No_HP_GFlow.xml

    # create a FIFO, write to it in the background
    PIPE=./cfg_${CFG}.fifo
    rm -f "$PIPE"
    mkfifo "$PIPE"
    envsubst '${CFG}' < "$TEMPLATE" > "$PIPE" &

    BIN=/pscratch/sd/h/hari_8/nEDM_project_LANL/build_script_28Feb2025/install/chroma_lanl/bin/chroma
   # srun -N 8 -n 32 -c 16 --ntasks-per-node=4 $BIN -geom 1 2 2 8 -i "$PIPE" -o /pscratch/sd/h/hari_8/nEDM_project_LANL/output_3pt2pt_flow/out_qtop_ggg_allFLOW_cfg_${CFG}.xml > log2pt3pt_${CFG}
    srun -N 8 -n 32 -c 16 --ntasks-per-node=4 $BIN -geom 1 2 2 8 -i "$PIPE"  > log2pt3pt_${CFG}
}

# List of CFG values
CFG_VALUES=(19010)

# Loop through each CFG value and call run_job
for CFG in "${CFG_VALUES[@]}"; do
    run_job $CFG
done

# (optional) wait for the background envsubst to finish
wait

# cleanup your fifo
rm -f ./cfg_${CFG}.fifo

