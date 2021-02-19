function fdr_imgWrapper(pairN, runN)

%% Wrapper around fdr_img function

% basics
alpha = 0.05;

% we assume we have the folders/files at ~/Desktop/testConfounds
folder = '/home/adamb/Desktop/testConfounds/results_Second/';
file1 = [folder, 'pair', num2str(pairN), '/coupling_run', num2str(runN), '_1.mat'];
saveF1 = [folder, 'pair', num2str(pairN), '/coupling_run', num2str(runN), '_1_FDR.mat'];
file2 = [folder, 'pair', num2str(pairN), '/coupling_run', num2str(runN), '_2.mat'];
saveF2 = [folder, 'pair', num2str(pairN), '/coupling_run', num2str(runN), '_2_FDR.mat'];

% call couplingFMRI for both
load(file1);
[bFDR, RsqFDR, pFFDR, pQ, lastK] = fdr_img(b, pF, Rsq, alpha);
save(saveF1, 'bFDR', 'RsqFDR', 'pFFDR', 'pQ', 'lastK', '-v7');

load(file2);
[bFDR, RsqFDR, pFFDR, pQ, lastK] = fdr_img(b, pF, Rsq, alpha);
save(saveF2, 'bFDR', 'RsqFDR', 'pFFDR', 'pQ', 'lastK', '-v7');

disp([char(10), 'Done with pair ', num2str(pairN), ', run ', num2str(runN)]);

return
