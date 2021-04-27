#!/bin/bash -l
# declare a name for this job to be sample_job
#PBS -N discovery_test_f00589z
# request the queue (enter the possible names, if omitted, default is the default)
# if more then 600 jobs use the largeq
#PBS -q testq
# request 1 core on 1 node
# ensure you reserve enough cores for the projected memory usage
# figuring 4G/core
#PBS -l nodes=1:ppn=1
# request 0 hours and 10 minutes of wall time
#PBS -l walltime=00:10:00
# mail is sent to you when the job begins and when it exits or aborts
# you can use all or some or none. If you don't want email leave this
# and the following (#PBS -M) out of the script.
#PBS -m bea
# specify your email address
#PBS -M f00589z@dartmouth.edu
## By default, PBS scripts execute in your home directory, not the
## directory from which they were submitted. The following line
## places you in the directory from which the job was submitted.
#cd $PBS_O_WORKDIR
# run the program
module load python
python discovery_python_test_nonInteractive.py