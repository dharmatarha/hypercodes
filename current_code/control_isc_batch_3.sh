#!/bin/bash

# Name of the job
#SBATCH --job-name=slurm_ISC_mem_per_cpu

# Number of compute nodes
#SBATCH --nodes=1

# Number of cores per node
#SBATCH --ntasks-per-node=16

# Number of CPUs per task
#SBATCH --cpus-per-task=1

# Request memory
#SBATCH --mem=256G

# Walltime (job duration)
#SBATCH --time=04:00:00

# Name of partition
#SBATCH --partition=standard

# Email address
#SBATCH --mail-user=f00589z@dartmouth.edu

# Email notifications (comma-separated options: BEGIN,END,FAIL)
# #SBATCH --mail-type=BEGIN,END,FAIL # remove for now

while getopts ":n:c:d:p:l:v:f:o:a:w:s:" opt; do
  case $opt in
    n) numJobs="$OPTARG"
    ;;
    c) circShift="$OPTARG"
    ;;
    d) debug="$OPTARG"
    ;;
    p) permutations="$OPTARG"
    ;;
    l) parallel="$OPTARG"
    ;;
    v) verbosity="$OPTARG"
    ;;
    f) fisherZtran="$OPTARG"
    ;;
    o) normalize="$OPTARG"
    ;;
    a) alpha="$OPTARG"
    ;;
    w) twoTailed="$OPTARG"
    ;;
    e) detrend="$OPTARG"
    ;;
    s) saveOutput="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

# set defaults
${numJobs:=16} # number of nodes
${circShift:=True} # use circle shifting
${debug:=False} # boolean, if True, truncates all time series to a small, manageable size (debugTRs x debugVox) to speed things up
${permutations:=1000} # number of iterations in permutation test
${parallel:=True} # boolean, if True, use the joblib parallelization function (speeds things up considerably)
${verbosity:=50} # joblib verbosity -- don't go over 50 lest ye want to print like a million outputs and slow everything down
${fisherZtran:=True} # boolean, if True, fisher z transform correlation coefficients
${normalize:=True} # boolean, if True, z-score time series before computing correlations
${alpha:=0.05} # alpha to use for permutation tests
${twoTailed:=True} # boolean, if True, uses a two-tailed permutation test, if False, uses a right-tailed test
${detrend:=False} # boolean, if True, uses detrended data
${saveOutput:=True} # boolean, if True, saves some ISC output variables in a .pkl file

# if parallel == false, then should reset numNodes to 1
printf "\n*** CONTROL TASK ISC PARAMETERS ***\n"
printf "Requested nodes: %s\n" "$numJobs"
printf "Circle shifting: %s\n" "$circShift"
printf "Debugging: %s\n" "$debug"
printf "Permutations: %s\n" "$permutations"
printf "Parallel computing: %s\n" "$parallel"
printf "Verbosity: %s\n" "$verbosity"
printf "Fisher's z transform null distributions: %s\n" "$fisherZtran"
printf "Standardize EPI time series before computing ISC: %s\n" "$normalize"
printf "Permutation test alpha: %s\n" "$alpha"
printf "Use a two-tailed permutation test: %s\n" "$twoTailed"
printf "Use detrended data: %s\n" "$detrend"
printf "Save ISC data: %s\n" "$saveOutput"

# run the program
ssh f00589z@discovery7 ". ~/.bashrc ; conda activate hypeScanCentral ; python3 /dartfs-hpc/rc/home/z/f00589z/hyperscanning/control_tasks/control_task_scripts/control_task_ISC.py -n $numJobs -c $circShift -d $debug -p $permutations -l $parallel -v $verbosity -f $fisherZtran -o $normalize -a $alpha -w $twoTailed -e $detrend -s $saveOutput"