function meanFDR_all(pairNumbers)

%% Function to collect all speaker-listener FDR corrected coeff images
% And to average the Rsq maps
%
% Also performs these separately for Joint and Ind conditions

if ~isvector(pairNumbers)
    error('Eh, pairNumbers, vector, now');
end

% get lookup variable
load('runLookup.mat');

% other basics
voxelN = 69880;
% RsqCutoff = 0.04;

% we assume we have the folders/files at ~/Desktop/testConfounds
folder = '/home/adamb/Desktop/testConfounds/results_Second/';

% preallocate
bAllJoint = zeros(length(pairNumbers)*2, 13, voxelN);
bAllInd = zeros(length(pairNumbers)*2, 13, voxelN);
RsqAllJoint = zeros(length(pairNumbers)*2, voxelN);
RsqAllInd = zeros(length(pairNumbers)*2, voxelN);

counterJoint = 0;
counterInd = 0;
for pairNo = pairNumbers
    
    for runN = 1:2
        
        for fileN = 1:2
        
            % load data
            fileFDR = [folder, 'pair', num2str(pairNo), '/coupling_run', num2str(runN), '_',num2str(fileN) , '_FDR.mat'];
            load(fileFDR);

            % if current run is a joint story
            if runLookup.joint(pairNo,2) == runN
                counterJoint = counterJoint + 1;
                bAllJoint(counterJoint, :, :) = bFDR;
                RsqAllJoint(counterJoint,:) = RsqFDR;

            elseif runLookup.ind(pairNo,2) == runN
                counterInd = counterInd + 1;
                bAllInd(counterInd, :, :) = bFDR;
                RsqAllInd(counterInd,:) = RsqFDR;     

            end
            
        end
        
    end
    
end


%% Get average maps

RsqMeanJoint = mean(RsqAllJoint,1);
RsqMeanInd = mean(RsqAllInd,1);
% cutoff


%% Save out results

save([folder, 'MeanMapsFDR.mat'], 'bAllJoint', 'bAllInd', 'RsqAllJoint', 'RsqAllInd', 'RsqMeanJoint', 'RsqMeanInd', '-v7');

return
    