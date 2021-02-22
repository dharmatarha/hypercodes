function [realFits, pseudoFits, parameters, p, FDR, Q, sigVox, time] = dbic_cbs_rollingISC(dbicSubs, cbsSubs, winTRs, stepSize, maxT, voxelCoords, groupPermuts, fitPermuts)

%% Function to run whole ISC process from start to finish without saving data -- will probably be very slow...
%
% !!! Second dataset, CBS + DBIC !!!
% 
% !!! For DRZEUSS !!!
%
% Goes through DBIC and CBS participants and calls couplingFMRI on pair 
% level data (both pseudo and real pairs)
%Inputs
%1) dbicSubs: 
%2) cbsSubs: 
%3) winTRs: window size [TRs]
%4) stepSize: step size [TRs]
%5) maxT: maximum number of TRs to use in either direction for lag model
%6) voxelCoords: voxel coordinates (1:69880 is whole brain -- need to
%modify to incorporate ROI maps)
%7) groupPermuts: # permutations to use in evaluating group level fit
%measures
%8) fitPermuts: (optional) # permutations to use in evaluating 
%pairwise/voxelwise model fits
%
%Outputs
%1) realFits: structure containing model fit statistics for real pairs
%2) pseudoFits: structure containing model fit statistics for pseudo pairs
%3) parameters: structure containing input parameters
%4) p: voxelwise p-values [voxel x window x condition] (conditions: 1=all,
%2=ind, 3=joint)
%5) FDR: FDR output [voxel x window x condition] (conditions: 1=all,
%2=ind, 3=joint)
%6) Q: q-values [voxel x window x condition] (conditions: 1=all,
%2=ind, 3=joint)
%7) sigVox: indices of voxels that show significant effects at q=alpha
%8) time: structure tracking duration of analysis

if nargin < 8
    stats = 0;
end

%% Basics

fileName = ['dbic' num2str(min(dbicSubs)) '-' num2str(max(dbicSubs)) '_cbs' num2str(min(cbsSubs)) '-' num2str(max(cbsSubs)) '_win' num2str(winTRs) '_step' num2str(stepSize) '_maxT' num2str(maxT)];

disp([char(10) '%%%%%%%%%%%%%%%%%'])
disp([char(10) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'])
disp([char(10) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'])
disp([char(10) 'starting analysis for... ' fileName])
disp([char(10) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'])
disp([char(10) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'])
disp([char(10) '%%%%%%%%%%%%%%%%%'])

% get lookup tables
load('/flash/wheatley/adamb/matlab/pairings2.mat');
load('/flash/wheatley/adamb/matlab/runLookup2.mat');

% range of subjects to go through
% dbicSubjects = 2:5;
% cbsSubjects = 2:9;

%get number of pseudo pairs
realPairN = length(intersect(dbicSubs, cbsSubs));
pseudoPairN = length(dbicSubs)*length(cbsSubs) - realPairN;

%keep track of how long things take
tic;
time.startDate = clock;
time.pairs = NaN(realPairN + pseudoPairN,3); %per pair duration -- column 1 'toc' output, column 2 gives approximate duration of pair ISC in minutes, column 3 shows real/pseudo status (1/0)

%rolling window parameters
totTRs = 615;
numSteps = ceil((totTRs - winTRs + 1) / stepSize); %total number of steps based on winTR and stepSize
lastTr = winTRs + (numSteps - 1) * stepSize; %last TR this approach will analyze
TRsLeftOut = totTRs - lastTr; %number of TRs that will be left out by the current approach
if TRsLeftOut > 0 
    warning(['The last ' num2str(TRsLeftOut) ' TRs will be left out due to the window and step size!'])
end

% preallocate 
pRows = pseudoPairN * 2;
rRows = realPairN * 2;
pseudoFits.ind.b = zeros(pRows, 2*maxT+1, length(voxelCoords), numSteps);
pseudoFits.ind.Rsq = zeros(pRows, length(voxelCoords), numSteps);
pseudoFits.ind.F = zeros(pRows, length(voxelCoords), numSteps);
pseudoFits.ind.pF = zeros(pRows, length(voxelCoords), numSteps);
pseudoFits.joint.b = zeros(pRows, 2*maxT+1, length(voxelCoords), numSteps);
pseudoFits.joint.Rsq = zeros(pRows, length(voxelCoords), numSteps);
pseudoFits.joint.F = zeros(pRows, length(voxelCoords), numSteps);
pseudoFits.joint.pF = zeros(pRows, length(voxelCoords), numSteps);
realFits.ind.b = zeros(rRows, 2*maxT+1, length(voxelCoords), numSteps);
realFits.ind.Rsq = zeros(rRows, length(voxelCoords), numSteps);
realFits.ind.F = zeros(rRows, length(voxelCoords), numSteps);
realFits.ind.pF = zeros(rRows, length(voxelCoords), numSteps);
realFits.joint.b = zeros(rRows, 2*maxT+1, length(voxelCoords), numSteps);
realFits.joint.Rsq = zeros(rRows, length(voxelCoords), numSteps);
realFits.joint.F = zeros(rRows, length(voxelCoords), numSteps);
realFits.joint.pF = zeros(rRows, length(voxelCoords), numSteps);

%preallocate map of which pairs correspond to which rows for the fit data
pseudoFits.pairMap.header = {'dbic','cbs','cbsSpeaker'};
pseudoFits.pairMap.data = NaN(pRows,length(pseudoFits.pairMap.header));
realFits.pairMap.header = {'dbic','cbs','cbsSpeaker'};
realFits.pairMap.data = NaN(rRows,length(realFits.pairMap.header));


%% Loop over DBIC paarticipants
pseudoCounter = 1;
realCounter = 1;
pairCounter = 1;
for dbicPairN = dbicSubs
    
    disp([char(10), 'Starting DBIC subject ', num2str(dbicPairN)]);
    
    % get dbicSub and runN
    dbicSub = pairings2.lookup{dbicPairN, 2}; % e.g. 'sid0000522'
    dbicRunInd = runLookup2.ind(dbicPairN,2);
    dbicRunJoint = runLookup2.joint(dbicPairN,2);
    
    % get independent condition files
    dbicFileSpeaker.ind = ['/flash/wheatley/adamb/hyperscanning_DBIC_ses2/sub-', dbicSub,...
    '_fmriprep/fmriprep/sub-', dbicSub, '/ses-pair0', num2str(dbicPairN), '/func/',...
    'sub-', dbicSub, '_ses-pair0', num2str(dbicPairN), '_task-storytelling', num2str(dbicRunInd), '_run-0',... 
    num2str(dbicRunInd), '_bold_space-MNI152NLin2009cAsym_preproc_nuisRegr_speaker.mat']; %***JDNote: TRs x voxels
    if ~exist(dbicFileSpeaker.ind, 'file')
        error([char(10), 'Cannot find dbicFile at ', dbicFileSpeaker.ind]);
    end
    dbicFileListener.ind = ['/flash/wheatley/adamb/hyperscanning_DBIC_ses2/sub-', dbicSub,...
    '_fmriprep/fmriprep/sub-', dbicSub, '/ses-pair0', num2str(dbicPairN), '/func/',...
    'sub-', dbicSub, '_ses-pair0', num2str(dbicPairN), '_task-storytelling', num2str(dbicRunInd), '_run-0',... 
    num2str(dbicRunInd), '_bold_space-MNI152NLin2009cAsym_preproc_nuisRegr_listener.mat'];
    if ~exist(dbicFileListener.ind, 'file')
        error([char(10), 'Cannot find dbicFile at ', dbicFileListener.ind]);
    end
    
    % get joint condition files
    dbicFileSpeaker.joint = ['/flash/wheatley/adamb/hyperscanning_DBIC_ses2/sub-', dbicSub,...
    '_fmriprep/fmriprep/sub-', dbicSub, '/ses-pair0', num2str(dbicPairN), '/func/',...
    'sub-', dbicSub, '_ses-pair0', num2str(dbicPairN), '_task-storytelling', num2str(dbicRunJoint), '_run-0',... 
    num2str(dbicRunJoint), '_bold_space-MNI152NLin2009cAsym_preproc_nuisRegr_speaker.mat']; %***JDNote: TRs x voxels
    if ~exist(dbicFileSpeaker.joint, 'file')
        error([char(10), 'Cannot find dbicFile at ', dbicFileSpeaker.joint]);
    end
    dbicFileListener.joint = ['/flash/wheatley/adamb/hyperscanning_DBIC_ses2/sub-', dbicSub,...
    '_fmriprep/fmriprep/sub-', dbicSub, '/ses-pair0', num2str(dbicPairN), '/func/',...
    'sub-', dbicSub, '_ses-pair0', num2str(dbicPairN), '_task-storytelling', num2str(dbicRunJoint), '_run-0',... 
    num2str(dbicRunJoint), '_bold_space-MNI152NLin2009cAsym_preproc_nuisRegr_listener.mat'];
    if ~exist(dbicFileListener.joint, 'file')
        error([char(10), 'Cannot find dbicFile at ', dbicFileListener.joint]);
    end
    
    % load independent condition data
    load(dbicFileSpeaker.ind); % results in variable dbicSpeaker
    load(dbicFileListener.ind); % results in variable dbicListener
    ind.dbicSpeaker = dbicSpeaker;
    ind.dbicListener = dbicListener;
    clear dbicSpeaker dbicListener
    
    % load joint condition data
    load(dbicFileSpeaker.joint); % results in variable dbicSpeaker
    load(dbicFileListener.joint); % results in variable dbicListener
    joint.dbicSpeaker = dbicSpeaker;
    joint.dbicListener = dbicListener;
    clear dbicSpeaker dbicListener
    
    disp([char(10), 'Loaded data for DBIC subject ', num2str(dbicPairN)]);
    
    
    %% Loop over CBS participants
    
    for cbsPairN = cbsSubs
        
        pairType = 0; %for adding to time.pairs below
         
        %%%%%%%%%%%%%%%%%
        %%% load data %%%
        %%%%%%%%%%%%%%%%%
        
        disp([char(10), 'Starting CBS subject ', num2str(cbsPairN)]);

        % get runN
        cbsRun = runLookup2.ind(cbsPairN,2);

        % get independent condition files
        cbsFileSpeaker.ind = ['/flash/wheatley/adamb/hyperscanning_CBS/sub-hid00000',...
            num2str(cbsPairN), '_fmriprep/fmriprep/sub-hid00000', num2str(cbsPairN),...
            '/ses-pair0', num2str(cbsPairN),'/func/sub-hid00000', num2str(cbsPairN), '_ses-pair0',...
            num2str(cbsPairN),'_task-storytelling', num2str(cbsRun), '_',...
            'run-0', num2str(cbsRun), '_bold_space-MNI152NLin2009cAsym_', ...
            'preproc_nuisRegr_speaker.mat'];
        if ~exist(cbsFileSpeaker.ind, 'file')
            error([char(10), 'Cannot find cbsFile at ', cbsFileSpeaker.ind]);
        end
        cbsFileListener.ind = ['/flash/wheatley/adamb/hyperscanning_CBS/sub-hid00000',...
            num2str(cbsPairN), '_fmriprep/fmriprep/sub-hid00000', num2str(cbsPairN),...
            '/ses-pair0', num2str(cbsPairN),'/func/sub-hid00000', num2str(cbsPairN), '_ses-pair0',...
            num2str(cbsPairN),'_task-storytelling', num2str(cbsRun), '_',...
            'run-0', num2str(cbsRun), '_bold_space-MNI152NLin2009cAsym_', ...
            'preproc_nuisRegr_listener.mat'];
        if ~exist(cbsFileListener.ind, 'file')
            error([char(10), 'Cannot find cbsFile at ', cbsFileListener.ind]);
        end

        % get joint condition files
        cbsFileSpeaker.joint = ['/flash/wheatley/adamb/hyperscanning_CBS/sub-hid00000',...
            num2str(cbsPairN), '_fmriprep/fmriprep/sub-hid00000', num2str(cbsPairN),...
            '/ses-pair0', num2str(cbsPairN),'/func/sub-hid00000', num2str(cbsPairN), '_ses-pair0',...
            num2str(cbsPairN),'_task-storytelling', num2str(cbsRun), '_',...
            'run-0', num2str(cbsRun), '_bold_space-MNI152NLin2009cAsym_', ...
            'preproc_nuisRegr_speaker.mat'];
        if ~exist(cbsFileSpeaker.joint, 'file')
            error([char(10), 'Cannot find cbsFile at ', cbsFileSpeaker.joint]);
        end
        cbsFileListener.joint = ['/flash/wheatley/adamb/hyperscanning_CBS/sub-hid00000',...
            num2str(cbsPairN), '_fmriprep/fmriprep/sub-hid00000', num2str(cbsPairN),...
            '/ses-pair0', num2str(cbsPairN),'/func/sub-hid00000', num2str(cbsPairN), '_ses-pair0',...
            num2str(cbsPairN),'_task-storytelling', num2str(cbsRun), '_',...
            'run-0', num2str(cbsRun), '_bold_space-MNI152NLin2009cAsym_', ...
            'preproc_nuisRegr_listener.mat'];
        if ~exist(cbsFileListener.joint, 'file')
            error([char(10), 'Cannot find cbsFile at ', cbsFileListener.joint]);
        end

        % load independent condition data
        load(cbsFileSpeaker.ind); % results in variable cbsSpeaker
        load(cbsFileListener.ind); % results in variable cbsListener
        ind.cbsSpeaker = cbsSpeaker;
        ind.cbsListener = cbsListener;
        clear cbsSpeaker cbsListener

        % load joint condition data
        load(cbsFileSpeaker.joint); % results in variable cbsSpeaker
        load(cbsFileListener.joint); % results in variable cbsListener
        joint.cbsSpeaker = cbsSpeaker;
        joint.cbsListener = cbsListener;
        clear cbsSpeaker cbsListener
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        %%% run couplingFMRI %%%
        %%%%%%%%%%%%%%%%%%%%%%%%
    
        if ~isequal(dbicPairN, cbsPairN) %if we're dealing with a pseudo pair...
            
            disp([char(10) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'])
            disp([char(10) '%%%%%%%%%%%%%% PSEUDO PAIR %%%%%%%%%%%%%%'])
            disp([char(10) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'])
            
            %% Loop over rolling windows...
            ROW = [(pseudoCounter*2-1) pseudoCounter*2];
            pseudoFits.pairMap.data(ROW,1) = dbicPairN;
            pseudoFits.pairMap.data(ROW,2) = cbsPairN;
            pseudoFits.pairMap.data(ROW,3) = [1; 0]; %first row has cbs as speaker
            for STEP = 1:numSteps %for each window...

                %get TRs to use in current window
                TRs = (1 + stepSize*(STEP-1)):(stepSize*(STEP-1) + winTRs);

                disp([char(10) 'Computing coupling fit measures for TRs ' num2str(1 + stepSize*(STEP-1)) ':' num2str(stepSize*(STEP-1) + winTRs)])

                %run couplingFMRI - independent condition
                [pseudoFits.ind.b(ROW(1),:,:,STEP), pseudoFits.ind.Rsq(ROW(1),:,STEP), pseudoFits.ind.F(ROW(1),:,STEP), pseudoFits.ind.pF(ROW(1),:,STEP)] = couplingFMRI(ind.cbsSpeaker(TRs,voxelCoords), ind.dbicListener(TRs,voxelCoords), maxT, fitPermuts);
                [pseudoFits.ind.b(ROW(2),:,:,STEP), pseudoFits.ind.Rsq(ROW(2),:,STEP), pseudoFits.ind.F(ROW(2),:,STEP), pseudoFits.ind.pF(ROW(2),:,STEP)] = couplingFMRI(ind.dbicSpeaker(TRs,voxelCoords), ind.cbsListener(TRs,voxelCoords), maxT, fitPermuts);
                
                %run couplingFMRI - joint condition
                [pseudoFits.joint.b(ROW(1),:,:,STEP), pseudoFits.joint.Rsq(ROW(1),:,STEP), pseudoFits.joint.F(ROW(1),:,STEP), pseudoFits.joint.pF(ROW(1),:,STEP)] = couplingFMRI(joint.cbsSpeaker(TRs,voxelCoords), joint.dbicListener(TRs,voxelCoords), maxT, fitPermuts);
                [pseudoFits.joint.b(ROW(2),:,:,STEP), pseudoFits.joint.Rsq(ROW(2),:,STEP), pseudoFits.joint.F(ROW(2),:,STEP), pseudoFits.joint.pF(ROW(2),:,STEP)] = couplingFMRI(joint.dbicSpeaker(TRs,voxelCoords), joint.cbsListener(TRs,voxelCoords), maxT, fitPermuts);

            end
            pseudoCounter = pseudoCounter + 1;
            
        
        else %if we're dealing with a real pair...
            
            pairType = 1; %for adding to time.pairs below
            
            disp([char(10) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'])
            disp([char(10) '%%%%%%%%%%%%%%% REAL PAIR %%%%%%%%%%%%%%%'])
            disp([char(10) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'])
            
            
            %% Loop over rolling windows...
            ROW = [(realCounter*2-1) realCounter*2];
            realFits.pairMap.data(ROW,1) = dbicPairN;
            realFits.pairMap.data(ROW,2) = cbsPairN;
            realFits.pairMap.data(ROW,3) = [1; 0]; %first row has cbs as speaker
            for STEP = 1:numSteps %for each window...

                %get TRs to use in current window
                TRs = (1 + stepSize*(STEP-1)):(stepSize*(STEP-1) + winTRs);

                disp([char(10) 'Computing coupling fit measures for TRs ' num2str(1 + stepSize*(STEP-1)) ':' num2str(stepSize*(STEP-1) + winTRs)])

                %run couplingFMRI - independent condition
                [realFits.ind.b(ROW(1),:,:,STEP), realFits.ind.Rsq(ROW(1),:,STEP), realFits.ind.F(ROW(1),:,STEP), realFits.ind.pF(ROW(1),:,STEP)] = couplingFMRI(ind.cbsSpeaker(TRs,voxelCoords), ind.dbicListener(TRs,voxelCoords), maxT, fitPermuts);
                [realFits.ind.b(ROW(2),:,:,STEP), realFits.ind.Rsq(ROW(2),:,STEP), realFits.ind.F(ROW(2),:,STEP), realFits.ind.pF(ROW(2),:,STEP)] = couplingFMRI(ind.dbicSpeaker(TRs,voxelCoords), ind.cbsListener(TRs,voxelCoords), maxT, fitPermuts);
                
                %run couplingFMRI - joint condition
                [realFits.joint.b(ROW(1),:,:,STEP), realFits.joint.Rsq(ROW(1),:,STEP), realFits.joint.F(ROW(1),:,STEP), realFits.joint.pF(ROW(1),:,STEP)] = couplingFMRI(joint.cbsSpeaker(TRs,voxelCoords), joint.dbicListener(TRs,voxelCoords), maxT, fitPermuts);
                [realFits.joint.b(ROW(2),:,:,STEP), realFits.joint.Rsq(ROW(2),:,STEP), realFits.joint.F(ROW(2),:,STEP), realFits.joint.pF(ROW(2),:,STEP)] = couplingFMRI(joint.dbicSpeaker(TRs,voxelCoords), joint.cbsListener(TRs,voxelCoords), maxT, fitPermuts);

            end
            realCounter = realCounter + 1;
            
            
            
        end
        disp([char(10), 'Done with CBS subject ', num2str(cbsPairN), ', got coupling measures']);
        
        %log time
        time.pairs(pairCounter,1) = toc;
        if pairCounter == 1
            time.pairs(pairCounter,2) = time.pairs(pairCounter,1);
        else
            time.pairs(pairCounter,2) = time.pairs(pairCounter,1) - time.pairs(pairCounter-1,1);
        end
        time.pairs(pairCounter,2) = pairType;
        pairCounter = pairCounter + 1;
        
        disp([char(10) '***************************************'])
        disp([char(10) '*** avg duration per pair: ' num2str(round(nanmean(time.pairs(:,1)) / 60,1))  ' min ***'])
        disp([char(10) '***************************************'])
        
        
         
    end
       
    disp([char(10), 'Done with DBIC subject ', num2str(dbicPairN)]);
        
    
end

disp([char(10), 'Model fitting finished!']);


%Permutation test to evaluate group level significance
disp([char(10) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'])
disp([char(10) '%%% Starting group level permutation tests %%%'])
disp([char(10) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'])

%Get some R^2 summary measures
for STEP = 1:numSteps
    
    %consolidate real pair data
    realFits.all.Rsq(:,:,STEP) = [realFits.joint.Rsq(:,:,STEP); realFits.ind.Rsq(:,:,STEP)];
    
    %consolidate pseudo pair data
    pseudoFits.all.Rsq(:,:,STEP) = [pseudoFits.joint.Rsq(:,:,STEP); pseudoFits.ind.Rsq(:,:,STEP)];
    
    %get mean across pairs
    realFits.all.Rsq_mean(STEP,:) = mean(realFits.all.Rsq(:,:,STEP), 1); %across independent and joint conditions
    realFits.ind.Rsq_mean(STEP,:) = mean(realFits.ind.Rsq(:,:,STEP), 1); %independent condition
    realFits.joint.Rsq_mean(STEP,:) = mean(realFits.joint.Rsq(:,:,STEP), 1); %joint condition
    
    %get SD across pairs
    realFits.all.Rsq_SD(STEP,:) = std(realFits.all.Rsq(:,:,STEP), 1); %across independent and joint conditions
    realFits.ind.Rsq_SD(STEP,:) = std(realFits.ind.Rsq(:,:,STEP), 1); %independent condition
    realFits.joint.Rsq_SD(STEP,:) = std(realFits.joint.Rsq(:,:,STEP), 1); %joint condition
    
    %get mean/SD across pairs then voxels
    realFits.all.Rsq_crossVox_mean(STEP,1) = mean(realFits.all.Rsq_mean(STEP,:)); %all
    realFits.all.Rsq_crossVox_SD(STEP,1) = std(realFits.all.Rsq_mean(STEP,:)); %all
    realFits.ind.Rsq_crossVox_mean(STEP,1) = mean(realFits.ind.Rsq_mean(STEP,:)); %ind
    realFits.ind.Rsq_crossVox_SD(STEP,1) = std(realFits.ind.Rsq_mean(STEP,:)); %ind
    realFits.joint.Rsq_crossVox_mean(STEP,1) = mean(realFits.joint.Rsq_mean(STEP,:)); %joint
    realFits.joint.Rsq_crossVox_SD(STEP,1) = std(realFits.joint.Rsq_mean(STEP,:)); %joint
    
end

%preallocate
p = NaN(numSteps,length(voxelCoords),3);

%% Loop through different tests (1=All, 2=Ind, 3=Joint)
for TEST = 1:3
    

    %% Loop through rolling windows...
    for STEP = 1:numSteps

        % test statistic is mean Rsq
        realDiff(STEP,:) = realFits.all.Rsq_mean(STEP,:,1) - mean(pseudoFits.all.Rsq(:,:,STEP), 1);

        % overall data
        if TEST == 1
            data = [realFits.all.Rsq(:,:,STEP); pseudoFits.all.Rsq(:,:,STEP)]; %get all pair data
            realN = size(realFits.all.Rsq,1); %number of real pairs
            pseudoN = size(pseudoFits.all.Rsq,1); %number of pseudo pairs
        elseif TEST == 2
            data = [realFits.ind.Rsq(:,:,STEP); pseudoFits.ind.Rsq(:,:,STEP)]; %get pair data from the independent condition
            realN = size(realFits.ind.Rsq,1); %number of real pairs
            pseudoN = size(pseudoFits.ind.Rsq,1); %number of pseudo pairs
        elseif TEST == 3
            data = [realFits.joint.Rsq(:,:,STEP); pseudoFits.joint.Rsq(:,:,STEP)]; %get pair data from the independent condition
            realN = size(realFits.joint.Rsq,1); %number of real pairs
            pseudoN = size(pseudoFits.joint.Rsq,1); %number of pseudo pairs
        end

        % preallocate results
        voxelN = size(realFits.all.Rsq,2);
        diffNull = zeros(groupPermuts, voxelN);

        % Loop
        for i = 1:groupPermuts

            %random mask
            mask = zeros(realN+pseudoN,1);
            mask(randperm(realN+pseudoN, realN)) = 1;

            %get null difference
            diffNull(i,:) = mean(data(mask==1,:),1)-mean(data(mask==0,:),1);

            %feedback
            if mod (i, 1000) == 0
                disp([char(10), 'Done with ', num2str(i)]);
            end

        end

        %% Get p values
        for v = 1:voxelN
            p(STEP,v,TEST) = (sum(diffNull(:,v) >= realDiff(v)) / groupPermuts);
        end
    end
end


%apply FDR correction and find voxels where a significant effect was observed (across tests/windows)
sigVox = [];
for TEST = 1:3
    for STEP = 1:numSteps
        [FDR(STEP,:,TEST),Q(STEP,:,TEST)] = mafdr(p(STEP,:,TEST)); %compare with Adam's FDR correction
        sigVox = [sigVox find(Q(STEP,:,TEST) <= 0.05)];
    end
end
sigVox = unique(sigVox);

%repeat plot above but only for significant voxels
%Get mean R^2 in each voxel across participants
for STEP = 1:numSteps
    
    %get mean across pairs
    realFits.all.Rsq_mean_SV(STEP,:) = mean(realFits.all.Rsq(:,sigVox,STEP), 1); %across independent and joint conditions
    realFits.ind.Rsq_mean_SV(STEP,:) = mean(realFits.ind.Rsq(:,sigVox,STEP), 1); %independent condition
    realFits.joint.Rsq_mean_SV(STEP,:) = mean(realFits.joint.Rsq(:,sigVox,STEP), 1); %joint condition
    
    %get SD across pairs
    realFits.all.Rsq_SD_SV(STEP,:) = std(realFits.all.Rsq(:,sigVox,STEP), 1); %across independent and joint conditions
    realFits.ind.Rsq_SD_SV(STEP,:) = std(realFits.ind.Rsq(:,sigVox,STEP), 1); %independent condition
    realFits.joint.Rsq_SD_SV(STEP,:) = std(realFits.joint.Rsq(:,sigVox,STEP), 1); %joint condition
    
    %get mean/SD across pairs then voxels
    realFits.all.Rsq_crossVox_mean_SV(STEP,1) = mean(realFits.all.Rsq_mean_SV(STEP,:)); %all
    realFits.all.Rsq_crossVox_SD_SV(STEP,1) = std(realFits.all.Rsq_mean_SV(STEP,:)); %all
    realFits.ind.Rsq_crossVox_mean_SV(STEP,1) = mean(realFits.ind.Rsq_mean_SV(STEP,:)); %ind
    realFits.ind.Rsq_crossVox_SD_SV(STEP,1) = std(realFits.ind.Rsq_mean_SV(STEP,:)); %ind
    realFits.joint.Rsq_crossVox_mean_SV(STEP,1) = mean(realFits.joint.Rsq_mean_SV(STEP,:)); %joint
    realFits.joint.Rsq_crossVox_SD_SV(STEP,1) = std(realFits.joint.Rsq_mean_SV(STEP,:)); %joint
    
end

yLimits = [0,0.3];

colors = {'k','r','b'};
figure('color','white'); 
subplot(1,2,1); hold on;
jitter = [-numSteps/50,0,numSteps/50];

h(1) = errorbar((1:numSteps) + jitter(TEST),realFits.all.Rsq_crossVox_mean',realFits.all.Rsq_crossVox_SD',['-o',colors{1}],'linewidth',2);
h(2) = errorbar((1:numSteps) + jitter(TEST),realFits.ind.Rsq_crossVox_mean',realFits.ind.Rsq_crossVox_SD',['-o',colors{2}],'linewidth',2);
h(3) = errorbar((1:numSteps) + jitter(TEST),realFits.joint.Rsq_crossVox_mean',realFits.joint.Rsq_crossVox_SD',['-o',colors{3}],'linewidth',2);

set(gca,'XTick', 1:numSteps,'FontSize',20);
ylim(yLimits)
xlabel('window','FontSize',20);
ylabel('mean R^2 (across subjects and voxels)','FontSize',20)
legend([h(1);h(2);h(3)],{'all (mean +/- SD across subs & voxels)','ind','joint'})


subplot(1,2,2); hold on;
colors = {'k','r','b'};
jitter = [-numSteps/50,0,numSteps/50];

h(1) = errorbar((1:numSteps) + jitter(TEST),realFits.all.Rsq_crossVox_mean_SV',realFits.all.Rsq_crossVox_SD_SV',['-o',colors{1}],'linewidth',2);
h(2) = errorbar((1:numSteps) + jitter(TEST),realFits.ind.Rsq_crossVox_mean_SV',realFits.ind.Rsq_crossVox_SD_SV',['-o',colors{2}],'linewidth',2);
h(3) = errorbar((1:numSteps) + jitter(TEST),realFits.joint.Rsq_crossVox_mean_SV',realFits.joint.Rsq_crossVox_SD_SV',['-o',colors{3}],'linewidth',2);

set(gca,'XTick', 1:numSteps,'FontSize',20);
ylim(yLimits)
xlabel('window','FontSize',20);
ylabel('mean R^2 (across subjects and voxels)','FontSize',20)
title('sig voxels')

%put parameters into structure for saving
parameters.dbicSubs = dbicSubs;
parameters.cbsSubs = cbsSubs;
parameters.winTRs = winTRs;
parameters.stepSize = stepSize;
parameters.maxT = maxT;
parameters.voxelCoords = voxelCoords;
parameters.groupPermuts = groupPermuts;
parameters.fitPermuts = fitPermuts;

time.totalDur = toc;
time.endData = clock;
disp([char(10) 'results file: /afs/.dbic.dartmouth.edu/usr/wheatley/jd/matlab/CBS_DBIC_Rolling/' fileName])
save(['/afs/.dbic.dartmouth.edu/usr/wheatley/jd/matlab/CBS_DBIC_Rolling/' fileName], 'realFits', 'parameters', 'p','FDR','Q','sigVox','time');

return
    
    
    
    
    
    
    
    