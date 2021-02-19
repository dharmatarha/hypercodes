function meanRsqImg_First(permut)

%% Function to calculate signifance of mean Rsq maps for First dataset
%
%
% !!! First dataset, DHMC + DBIC !!!!
%
% !!! DRZEUSS !!!
%
%
%

if nargin ~= 1
    error('Need input arg "permut", the number of premutations to perform');
end


%% Basics

voxelN = 69880;


%% Load real pair data

folder1 = '/flash/wheatley/adamb/matlab/DHMC_DBIC/';
load([folder1, 'AllPairData.mat']);

% relevant matrices
RsqAll = [RsqAllJoint;RsqAllInd];

disp([char(10), 'Got real data, with size ', num2str(size(RsqAll,1)), ' X ', num2str(size(RsqAll,2))]);


%% Load pseudo data

folder2 = '/flash/wheatley/adamb/matlab/permFirst/';
load([folder2, 'permFirst_Joint_1.mat']);
Rsq1 = Rsq; clearvars Rsq;
load([folder2, 'permFirst_Joint_2.mat']);
Rsq2 = Rsq; clearvars Rsq;
load([folder2, 'permFirst_Ind_1.mat']);
Rsq3 = Rsq; clearvars Rsq;
load([folder2, 'permFirst_Ind_2.mat']);
Rsq4 = Rsq; clearvars Rsq;

RsqJointPseudo = [Rsq1; Rsq2];
RsqIndPseudo = [Rsq3; Rsq4];
RsqAllPseudo = [Rsq1; Rsq2; Rsq3; Rsq4];
clearvars Rsq1 Rsq2 Rsq3 Rsq4;

% save path
saveFile = [folder1, 'permResults.mat'];

disp([char(10), 'Got pseudo data, with size ', num2str(size(RsqAllPseudo,1)),...
    ' X ', num2str(size(RsqAllPseudo,2))]);


%% Permutation test for grand average

% number of real and pseudo data couples
realN = size(RsqAll,1);
pseudoN = size(RsqAllPseudo,1);

% test statistic is mean Rsq
realDiff = mean(RsqAll, 1)-mean(RsqAllPseudo,1);

% overall data
data = [RsqAll;RsqAllPseudo];

% preallocate results
diffNull = zeros(permut, voxelN);

% Loop
for i = 1:permut
    
    % random mask
    mask = zeros(realN+pseudoN,1);
    mask(randperm(realN+pseudoN, realN)) = 1;
    
    % get test stat for the two groups
%     group = mean(data(mask==1,:),1);
%     rest = mean(data(mask==0,:),1);
    
    diffNull(i,:) = mean(data(mask==1,:),1)-mean(data(mask==0,:),1);
    
    % feedback
    
    if mod (i, 1000) == 0
        disp([char(10), 'Done with ', num2str(i)]);
    end
    
    
end    
    
RsqImg = mean(RsqAll, 1);


%% Get p values
p = zeros(1, voxelN);

for v = 1:voxelN
    
    p(v) = (sum(diffNull(:,v) >= realDiff(v)) / permut);
    
end
    

save(saveFile, 'RsqImg', 'p');

return
        
    
    
    
    


