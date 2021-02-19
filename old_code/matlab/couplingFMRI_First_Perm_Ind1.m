function couplingFMRI_First_Perm_Ind1

%% Function to create pseudo-pair coupling data
%
% !!! First dataset, DHMC + DBIC !!!
% 
% !!! For DRZEUSS !!!
%
% Goes through DBIC participants 1-5 (2:6), pairs them up with all
% potential DHMC participants (except the real partner) and call
% couplingFMRI on their pair level data.
%
% Only Ind conditions!
%
% Saves out results to /flash.../matlab/permFirst/...
%


%% Basics

% get lookup tables
load('/flash/wheatley/adamb/matlab/pairings.mat');
load('/flash/wheatley/adamb/matlab/runLookup.mat');

% range of subjects to go through
dbicSubjects = 2:6;
dhmcSubjects = [2:8, 10:12];
pseudoPairN = length(dbicSubjects)*length(dhmcSubjects);

% params for couplingFMRI
maxT = 3;
stats = [];

% misc params
voxelN = 69880;

% preallocate 
b = zeros(pseudoPairN, 2*maxT+1, voxelN);
Rsq = zeros(pseudoPairN, voxelN);
F = zeros(pseudoPairN, voxelN);
pF = zeros(pseudoPairN, voxelN);

% save file location
savePath = '/flash/wheatley/adamb/matlab/permFirst/permFirst_Ind_1.mat';


%% Loop over DBIC paarticipants
counter = 0;
for dbicPairN = dbicSubjects
    
    disp([char(10), 'Starting DBIC subject ', num2str(dbicPairN)]);
    
    % get dbicSub and runN
    dbicSub = pairings.lookup{dbicPairN, 2}; % e.g. 'sid0000522'
    dbicRun = runLookup.ind(dbicPairN,2);

    % get files
    dbicFileSpeaker = ['/flash/wheatley/adamb/hyperscanning_DBIC/sub-', dbicSub,...
    '_fmriprep/fmriprep/sub-', dbicSub, '/func/',...
    'sub-', dbicSub, '_task-storytelling', num2str(dbicRun), '_run-0',... 
    num2str(dbicRun), '_bold_space-MNI152NLin2009cAsym_preproc_nuisRegr_speaker.mat'];
    if ~exist(dbicFileSpeaker, 'file')
        error([char(10), 'Cannot find dbicFileSpeaker at ', dbicFileSpeaker]);
    end
    dbicFileListener = ['/flash/wheatley/adamb/hyperscanning_DBIC/sub-', dbicSub,...
    '_fmriprep/fmriprep/sub-', dbicSub, '/func/',...
    'sub-', dbicSub, '_task-storytelling', num2str(dbicRun), '_run-0',... 
    num2str(dbicRun), '_bold_space-MNI152NLin2009cAsym_preproc_nuisRegr_listener.mat'];
    if ~exist(dbicFileListener, 'file')
        error([char(10), 'Cannot find dbicFileSpeaker at ', dbicFileListener]);
    end
    
    % load data
    load(dbicFileSpeaker); % results in variable dbicSpeaker
    load(dbicFileListener); % results in variable dbicListener
    
    disp([char(10), 'Loaded data for DBIC subject ', num2str(dbicPairN)]);
    
    
    %% Loop over DHMC participants
    
    for dhmcPairN = dhmcSubjects
        
        % only go on if pseudo pair
        if ~isequal(dbicPairN, dhmcPairN)
            
            disp([char(10), 'Starting DHMC subject ', num2str(dhmcPairN)]);
            
            % get runN
            dhmcRun = runLookup.ind(dhmcPairN,2);

            % get files
            dhmcFileSpeaker = ['/flash/wheatley/adamb/hyperscanning_DHMC/sub-pair',...
                num2str(dhmcPairN), 'DHMC_fmriprep/fmriprep/sub-pair', num2str(dhmcPairN),...
                'DHMC/func/sub-pair', num2str(dhmcPairN), 'DHMC_task-storytelling_',...
                'acq-3mm_run-', num2str(dhmcRun), '_bold_space-MNI152NLin2009cAsym_', ...
                'preproc_nuisRegr_speaker.mat'];
            if ~exist(dhmcFileSpeaker, 'file')
                error([char(10), 'Cannot find dhmcFile at ', dhmcFileSpeaker]);
            end
            dhmcFileListener = ['/flash/wheatley/adamb/hyperscanning_DHMC/sub-pair',...
                num2str(dhmcPairN), 'DHMC_fmriprep/fmriprep/sub-pair', num2str(dhmcPairN),...
                'DHMC/func/sub-pair', num2str(dhmcPairN), 'DHMC_task-storytelling_',...
                'acq-3mm_run-', num2str( dhmcRun), '_bold_space-MNI152NLin2009cAsym_', ...
                'preproc_nuisRegr_listener.mat'];
            if ~exist(dhmcFileListener, 'file')
                error([char(10), 'Cannot find dhmcFile at ', dhmcFileListener]);
            end

            % load data
            load(dhmcFileSpeaker); % results in variable dhmcSpeaker
            load(dhmcFileListener); % results in variable dhmcListener
    
            % call couplingFMRI, no stats param
            counter = counter+1;
            [b(counter,:,:), Rsq(counter,:), F(counter,:), pF(counter,:)] = couplingFMRI(dhmcSpeaker, dbicListener, maxT);
            counter = counter+1;
            [b(counter,:,:), Rsq(counter,:), F(counter,:), pF(counter,:)] = couplingFMRI(dbicSpeaker, dhmcListener, maxT);
            
            disp([char(10), 'Donoe with DHMC subject ', num2str(dhmcPairN), ', got coupling measures']);
            
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
    
    
    
    
    
    
    
    