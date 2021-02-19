function meanNonFDR_Second_local(pairNumbers)

%% Function to collect all speaker-listener FDR corrected coeff images
% Collect all coupling outputs into one mat file
%
% !!! Second Dataset, DRZEUSS !!!
%
%


if ~isvector(pairNumbers)
    error('Eh, pairNumbers, vector, now');
end

% get lookup variable
load('runLookup2.mat');

% other basics
voxelN = 69880;

% we assume we have the files under /flash..../matlab/CBS_DBIC/pairX/
folder = '/home/adamb/Desktop/testConfounds/results_Second/';

% preallocate
bAllJoint = zeros(length(pairNumbers)*2, 13, voxelN);
bAllInd = zeros(length(pairNumbers)*2, 13, voxelN);
RsqAllJoint = zeros(length(pairNumbers)*2, voxelN);
RsqAllInd = zeros(length(pairNumbers)*2, voxelN);
RsqPairContrast = zeros(length(pairNumbers)*2, voxelN);

counterJoint = 0;
counterInd = 0;
for pairNo = pairNumbers
    
    for runN = 1:2
        
        for fileN = 1:2
        
            % load data
            file = [folder, 'pair', num2str(pairNo), '/coupling_run', num2str(runN), '_',num2str(fileN) , '.mat'];
            load(file);

            % if current run is a joint story
            if runLookup2.joint(pairNo,2) == runN
                counterJoint = counterJoint + 1;
                bAllJoint(counterJoint, :, :) = b;
                RsqAllJoint(counterJoint,:) = Rsq;

            elseif runLookup2.ind(pairNo,2) == runN
                counterInd = counterInd + 1;
                bAllInd(counterInd, :, :) = b;
                RsqAllInd(counterInd,:) = Rsq;     

            end
            
        end
        
    end
    
end


%% Get average maps

RsqMeanJoint = mean(RsqAllJoint,1);
RsqMeanInd = mean(RsqAllInd,1);
RsqMeanAll = mean([RsqAllJoint;RsqAllInd],1);


%% Get pair-level contrasts
for i = 1:length(pairNumbers)*2
    RsqPairContrast(i,:) = RsqAllJoint(i,:) - RsqAllInd(i,:);
end
RsqMeanContrast = mean(RsqPairContrast, 1);


%% Save out results

save([folder, 'AllPairData_local.mat'], 'bAllJoint', 'bAllInd', 'RsqAllJoint',...
    'RsqAllInd', 'RsqMeanJoint', 'RsqMeanInd', 'RsqMeanAll', 'RsqPairContrast', '-v7');

return