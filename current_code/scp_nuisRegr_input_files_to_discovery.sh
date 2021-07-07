#!/bin/sh

# This script transfers the necessary input files for running nuisance regression on the hyperscanning control tasks
# (listening and reading) from the drzeuss server to the discovery server.

# define DBIC and CBS subject ID arrays
dbicIDs=("sid000007" "sid000009" "sid000560" "sid000535" "sid000102" "sid000416" "sid000499" "sid000142");
cbsIDs=("hid000002" "hid000003" "hid000004" "hid000005" "hid000006" "hid000007" "hid000008" "hid000009");

# transfer EPI .nii.gz and confound .csv files -- see Hyperscanning Hub for info on how these were generated:
# https://docs.google.com/document/d/1uptaoaXoeYfHznZ8dFk6LxUxm2cncoLP015GBpot2MU/edit#
for PAIR in {2..9}
do
  for RUN in {3..4}
  do
    # DBIC preproc.nii.gz
    fileName=sub-${dbicIDs[${PAIR}-2]}_ses-pair0${PAIR}_task-storytelling${RUN}_run-0${RUN}_bold_space-MNI152NLin2009cAsym_preproc.nii.gz
    rsync -av jd@drzeuss.dartmouth.edu:/flash/wheatley/adamb/hyperscanning_DBIC_ses2/sub-${dbicIDs[${PAIR}-2]}_fmriprep/fmriprep/sub-${dbicIDs[${PAIR}-2]}/ses-pair0${PAIR}/func/${fileName} /dartfs-hpc/rc/home/z/f00589z/hyperscanning/control_tasks/nuisRegr_input_files/${fileName}

    # DBIC confounds_truncated_2021.csv
    fileName=sub-${dbicIDs[${PAIR}-2]}_ses-pair0${PAIR}_task-storytelling${RUN}_run-0${RUN}_bold_confounds_truncated_2021.csv
    rsync -av jd@drzeuss.dartmouth.edu:/afs/dbic.dartmouth.edu/usr/wheatley/jd/control_tasks/${fileName} /dartfs-hpc/rc/home/z/f00589z/hyperscanning/control_tasks/nuisRegr_input_files/${fileName}

    # CBS preproc.nii.gz
    fileName=sub-${cbsIDs[${PAIR}-2]}_ses-pair0${PAIR}_task-storytelling${RUN}_run-0${RUN}_bold_space-MNI152NLin2009cAsym_preproc.nii.gz
    rsync -av jd@drzeuss.dartmouth.edu:/flash/wheatley/adamb/hyperscanning_CBS/sub-${cbsIDs[${PAIR}-2]}_fmriprep/fmriprep/sub-${cbsIDs[${PAIR}-2]}/ses-pair0${PAIR}/func/${fileName} /dartfs-hpc/rc/home/z/f00589z/hyperscanning/control_tasks/nuisRegr_input_files/${fileName}

    # CBS confounds_truncated_2021.csv
    fileName=sub-${cbsIDs[${PAIR}-2]}_ses-pair0${PAIR}_task-storytelling${RUN}_run-0${RUN}_bold_confounds_truncated_2021.csv
    rsync -av jd@drzeuss.dartmouth.edu:/afs/dbic.dartmouth.edu/usr/wheatley/jd/control_tasks/${fileName} /dartfs-hpc/rc/home/z/f00589z/hyperscanning/control_tasks/nuisRegr_input_files/${fileName}
  done
done

# copy over brain mask as well
rsync -av jd@drzeuss.dartmouth.edu:/flash/wheatley/adamb/mni_asym09c_mask_resamp3x3.nii.gz /dartfs-hpc/rc/home/z/f00589z/hyperscanning/control_tasks/nuisRegr_input_files/mni_asym09c_mask_resamp3x3.nii.gz