function couplingFMRI_Second_Perm_Ind2

%% Function to create pseudo-pair coupling data
%
% !!! Second dataset, CBS + DBIC !!!
% 
% !!! For DRZEUSS !!!
%
% Goes through DBIC participants 5-8 (6:9), pairs them up with all
% potential CBS participants (except the real partner) and call
% couplingFMRI on their pair level data.
%
% Only Ind conditions!
%
% Saves out results to /flash.../matlab/permSecond/...
%


%% Basics

% get lookup tables
load('/flash/wheatley/adamb/matlab/pairings2.mat');
load('/flash/wheatley/adamb/matlab/runLookup2.mat');

% range of subjects to go through
dbicSubjects = 6:9;
cbsSubjects = 2:9;
pseudoPairN = length(dbicSubjects)*(length(cbsSubjects)-1);

% params for couplingFMRI
maxT = 6;
stats = [];

% misc params
voxelN = 69880;

% preallocate 
b = zeros(pseudoPairN*2, 2*maxT+1, voxelN);
Rsq = zeros(pseudoPairN*2, voxelN);
F = zeros(pseudoPairN*2, voxelN);
pF = zeros(pseudoPairN*2, voxelN);

% save file location
savePath = '/flash/wheatley/adamb/matlab/permSecond/permSecond_Ind_2.mat';


%% Loop over DBIC paarticipants
counter = 0;
for dbicPairN = dbicSubjects
    
    disp([char(10), 'Starting DBIC subject ', num2str(dbicPairN)]);
    
    % get dbicSub and runN
    dbicSub = pairings2.lookup{dbicPairN, 2}; % e.g. 'sid0000522'
    dbicRun = runLookup2.ind(dbicPairN,2);

    % get files
    dbicFileSpeaker = ['/flash/wheatley/adamb/hyperscanning_DBIC_ses2/sub-', dbicSub,...
    '_fmriprep/fmriprep/sub-', dbicSub, '/ses-pair0', num2str(dbicPairN), '/func/',...
    'sub-', dbicSub, '_ses-pair0', num2str(dbicPairN), '_task-storytelling', num2str(dbicRun), '_run-0',... 
    num2str(dbicRun), '_bold_space-MNI152NLin2009cAsym_preproc_nuisRegr_speaker.mat'];
    if ~exist(dbicFileSpeaker, 'file')
        error([char(10), 'Cannot find dbicFileSpeaker at ', dbicFileSpeaker]);
    end
    dbicFileListener = ['/flash/wheatley/adamb/hyperscanning_DBIC_ses2/sub-', dbicSub,...
    '_fmriprep/fmriprep/sub-', dbicSub, '/ses-pair0', num2str(dbicPairN), '/func/',...
    'sub-', dbicSub, '_ses-pair0', num2str(dbicPairN), '_task-storytelling', num2str(dbicRun), '_run-0',... 
    num2str(dbicRun), '_bold_space-MNI152NLin2009cAsym_preproc_nuisRegr_listener.mat'];
    if ~exist(dbicFileListener, 'file')
        error([char(10), 'Cannot find dbicFileSpeaker at ', dbicFileListener]);
    end
    
    % load data
    load(dbicFileSpeaker); % results in variable dbicSpeaker
    load(dbicFileListener); % results in variable dbicListener
    
    disp([char(10), 'Loaded data for DBIC subject ', num2str(dbicPairN)]);
    
    
    %% Loop over CBS participants
    
    for cbsPairN = cbsSubjects
        
        % only go on if pseudo pair
        if ~isequal(dbicPairN, cbsPairN)
            
            disp([char(10), 'Starting CBS subject ', num2str(cbsPairN)]);
            
            % get runN
            cbsRun = runLookup2.ind(cbsPairN,2);

            % get files
            cbsFileSpeaker = ['/flash/wheatley/adamb/hyperscanning_CBS/sub-hid00000',...
                num2str(cbsPairN), '_fmriprep/fmriprep/sub-hid00000', num2str(cbsPairN),...
                '/ses-pair0', num2str(cbsPairN),'/func/sub-hid00000', num2str(cbsPairN), '_ses-pair0',...
                num2str(cbsPairN),'_task-storytelling', num2str(cbsRun), '_',...
                'run-0', num2str(cbsRun), '_bold_space-MNI152NLin2009cAsym_', ...
                'preproc_nuisRegr_speaker.mat'];
            if ~exist(cbsFileSpeaker, 'file')
                error([char(10), 'Cannot find cbsFile at ', cbsFileSpeaker]);
            end
            cbsFileListener = ['/flash/wheatley/adamb/hyperscanning_CBS/sub-hid00000',...
                num2str(cbsPairN), '_fmriprep/fmriprep/sub-hid00000', num2str(cbsPairN),...
                '/ses-pair0', num2str(cbsPairN),'/func/sub-hid00000', num2str(cbsPairN), '_ses-pair0',...
                num2str(cbsPairN),'_task-storytelling', num2str(cbsRun), '_',...
                'run-0', num2str(cbsRun), '_bold_space-MNI152NLin2009cAsym_', ...
                'preproc_nuisRegr_listener.mat'];
            if ~exist(cbsFileListener, 'file')
                error([char(10), 'Cannot find cbsFile at ', cbsFileListener]);
            end

            % load data
            load(cbsFileSpeaker); % results in variable cbsSpeaker
            load(cbsFileListener); % results in variable cbsListener
    
            % call couplingFMRI, no stats param
            counter = counter+1;
            [b(counter,:,:), Rsq(counter,:), F(counter,:), pF(counter,:)] = couplingFMRI(cbsSpeaker, dbicListener, maxT);
            counter = counter+1;
            [b(counter,:,:), Rsq(counter,:), F(counter,:), pF(counter,:)] = couplingFMRI(dbicSpeaker, cbsListener, maxT);
            
            disp([char(10), 'Done with CBS subject ', num2str(cbsPairN), ', got coupling measures']);
            
        end
        
    end
       
    
    % save out the end of every dbic subject
    save(savePath, 'b', 'Rsq', 'F', 'pF');
    
    disp([char(10), 'Done with DBIC subject ', num2str(dbicPairN)]);
        
    
end
    

% save again
save(savePath, 'b', 'Rsq', 'F', 'pF');

disp([char(10), 'Everything finished!']);

return
    
    
    
    
    
    
    
    