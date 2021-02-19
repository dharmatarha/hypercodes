function readConfoundsDHMC(subject)



results=tdfread(['sub-pair',num2str(subject),'DHMC_task-storytelling_acq-3mm_run-1_bold_confounds.tsv'],'\t');
save(['subject',num2str(subject),'DHMC_1.mat'],'results');
results=tdfread(['sub-pair',num2str(subject),'DHMC_task-storytelling_acq-3mm_run-2_bold_confounds.tsv'],'\t');
save(['subject',num2str(subject),'DHMC_2.mat'],'results');
results=tdfread(['sub-pair',num2str(subject),'DHMC_task-listening_acq-3mm_run-3_bold_confounds.tsv'],'\t');
save(['subject',num2str(subject),'DHMC_3.mat'],'results');
results=tdfread(['sub-pair',num2str(subject),'DHMC_task-reading_acq-3mm_run-4_bold_confounds.tsv'],'\t');
save(['subject',num2str(subject),'DHMC_4.mat'],'results');

