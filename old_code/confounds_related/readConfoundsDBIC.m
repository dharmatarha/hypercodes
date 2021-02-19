function readConfoundsDBIC(sid)



results=tdfread(['sub-sid000',sid,'_task-storytelling1_run-01_bold_confounds.tsv'],'\t');
save(['subject',sid,'_1.mat'],'results');
results=tdfread(['sub-sid000',sid,'_task-storytelling2_run-02_bold_confounds.tsv'],'\t');
save(['subject',sid,'_2.mat'],'results');
results=tdfread(['sub-sid000',sid,'_task-storytelling3_run-03_bold_confounds.tsv'],'\t');
save(['subject',sid,'_3.mat'],'results');
results=tdfread(['sub-sid000',sid,'_task-storytelling4_run-04_bold_confounds.tsv'],'\t');
save(['subject',sid,'_4.mat'],'results');

