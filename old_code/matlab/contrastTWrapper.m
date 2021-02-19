function contrastTWrapper

%% Quick and dirty wrapper function for contrastT
%
%


%% Basics, load data we use

contrastEarly = [0.5, 0.5, -0.2, -0.2, -0.2, -0.2, -0.2];
contrastSync = [-0.25, -0.25, 0.33, 0.33, 0.33, -0.25, -0.25];
contrastLate = [-0.2, -0.2, -0.2, -0.2, -0.2, 0.5, 0.5];

load('/home/adamb/Desktop/testConfounds/results_First/AllPairData.mat');
load('/home/adamb/Desktop/testConfounds/results_First/FDRImages.mat');

voxelMask=RsqMeanAllFDR>0;

bAll = [bAllJoint; bAllInd];

voxelN = length(voxelMask);

%% Call contrastT

[lEarly, pEarly, bEarly] = contrastT(bAll, voxelMask, contrastEarly);
[lSync, pSync, bSync] = contrastT(bAll, voxelMask, contrastSync);
[lLate, pLate, bLate] = contrastT(bAll, voxelMask, contrastLate);


%% Get FDR for each
[pQEarly, ~, ~, ~] = fdr_quick(bEarly, pEarly);
pEarly2 = pEarly; pEarly2(pEarly2>pQEarly)=nan;
disp([char(10), num2str(sum(~isnan(pEarly2))), ' voxels after FDR for pEarly']);

[pQSync, ~, ~, ~] = fdr_quick(bSync, pSync);
pSync2 = pSync; pSync2(pSync2>pQSync)=nan;
disp([char(10), num2str(sum(~isnan(pSync2))), ' voxels after FDR for pSync']);

[pQLate, ~, ~, ~] = fdr_quick(bLate, pLate);
pLate2 = pLate; pLate2(pLate2>pQLate)=nan;
disp([char(10), num2str(sum(~isnan(pLate2))), ' voxels after FDR for pLate']);


%% Get mask for image

mask = zeros(1, voxelN);
% detect overlaps as well - if there is any, take the value with smallest p
[Y, idx] = min([pEarly2', pSync2', pLate2'], [], 2, 'omitnan');
mask = idx; mask(isnan(Y))=nan; mask = mask';


%% Save things out

save('voxelLagsFDR.mat', 'lEarly', 'pEarly', 'pEarly2', 'bEarly', 'pQEarly', 'lSync', 'pSync', 'pSync2', 'bSync', 'pQSync', 'lLate', 'pLate', 'pLate2', 'bLate', 'pQLate', 'mask', '-v7');

return
