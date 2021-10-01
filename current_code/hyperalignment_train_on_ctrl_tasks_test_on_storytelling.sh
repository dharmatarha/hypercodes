#!/bin/bash

# Name of the job
#SBATCH --job-name=hypAlign_on_ctrl_tasks_test_on_storytelling

# Number of compute nodes
#SBATCH --nodes=1

# Number of cores per node
#SBATCH --ntasks-per-node=16

# Number of CPUs per task
#SBATCH --cpus-per-task=1

# Request memory
#SBATCH --mem=256G

# Walltime (job duration)
#SBATCH --time=12:00:00

# Name of partition
#SBATCH --partition=standard

# Email address
#SBATCH --mail-user=f00589z@dartmouth.edu

# Email notifications (comma-separated options: BEGIN,END,FAIL)
#SBATCH --mail-type=BEGIN,END,FAIL

source /optnfs/common/miniconda3/etc/profile.d/conda.sh
conda activate py27Env
python2.7 /dartfs-hpc/rc/lab/W/WheatleyT/f00589z/hyperscanning/preprocessing/hyperalignment/hyperalignment_train_on_ctrl_tasks_test_on_storytelling.py -r 5 -c 1 -e 0