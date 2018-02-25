#!/bin/bash
#
#SBATCH --job-name=mcl1t
#SBATCH --output=mcl1t.txt
#
#SBATCH --time=20:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=500
#SBATCH --gres gpu:8
#SBATCH --gres-flags=enforce-binding

#srun hostname
#srun sleep 600

echo "Number of nodes: $SLURM_NNODES"
echo "Nodelist: $SLURM_NODELIST"
echo "Number of tasks: $SLURM_NTASKS"
echo "Tasks per node: $SLURM_TASKS_PER_NODE"

scontrol show hostname $SLURM_NODELIST > .slurm_nodes
awk '{ for(i=0;i<8;i++)print $1 ",""0:"i",1,centos-OpenCL,,/tmp"}' < .slurm_nodes > nodefile

python ~/src/async_re-openmm/bedamtempt_async_re.py mcl1-23.cntl  > LOG 2>&1
