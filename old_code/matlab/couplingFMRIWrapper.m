function couplingFMRIWrapper(pairN, runN, stats)

%% Wrapper function around couplingFMRI.m
%
% !!! First dataset, DBIC + DHMC !!!
%
% Inputs are the usual pairN and runN
%
% Outputs are saved out files
%


%% Input checks

if nargin < 2
    error('Need input args pairN and runN');
end
if nargin < 3
    stats = [];
else
    if (stats < 0) || (stats > 10000)
        error([char(10), 'Input arg stats should be between 0 - 10000, eh?']);
    end
end
if ~ismember(pairN,1:14)
    error('Input arg pairN should be int between 1 and 14');
end
% we handle the special case of pair 14 as well
if ismember(pairN, 1:13)
    if ~ismember(runN,1:2)
        error('Input arg runN should be 1 or 2');
    end
elseif pairN == 14
    if ~ismember(runN,[1,2,5])
        error('Input arg runN should be 1 or 2 or 5'); % 5 is for repeated run1
    end
end

disp([char(10), 'Started couplingFMRIWrapper.m with inputs ',...
    num2str(pairN), ' (pair no) and ', num2str(runN), ' (run no)']);


%% Get file path to everything we need

% get the DBIC subject name for given pair
load('/flash/wheatley/adamb/matlab/pairings.mat');
dbicSub = pairings.lookup{pairN, 2}; % e.g. 'sid0000522'

% get path to DBIC epi files
dbicFileSpeaker = ['/flash/wheatley/adamb/hyperscanning_DBIC/sub-', dbicSub,...
    '_fmriprep/fmriprep/sub-', dbicSub, '/func/',...
    'sub-', dbicSub, '_task-storytelling', num2str(runN), '_run-0',... 
    num2str(runN), '_bold_space-MNI152NLin2009cAsym_preproc_nuisRegr_speaker.mat'];
if ~exist(dbicFileSpeaker, 'file')
    error([char(10), 'Cannot find dbicFileSpeaker at ', dbicFileSpeaker]);
end
% get parts of filename
[dbicFileSpeakerPath, dbicFileSpeakerBase, ~] = fileparts(dbicFileSpeaker);

dbicFileListener = ['/flash/wheatley/adamb/hyperscanning_DBIC/sub-', dbicSub,...
    '_fmriprep/fmriprep/sub-', dbicSub, '/func/',...
    'sub-', dbicSub, '_task-storytelling', num2str(runN), '_run-0',... 
    num2str(runN), '_bold_space-MNI152NLin2009cAsym_preproc_nuisRegr_listener.mat'];
if ~exist(dbicFileListener, 'file')
    error([char(10), 'Cannot find dbicFileSpeaker at ', dbicFileListener]);
end
% get parts of filename
[dbicFileListenerPath, dbicFileListenerBase, ~] = fileparts(dbicFileListener);


% get path to DHMC epi files
dhmcFileSpeaker = ['/flash/wheatley/adamb/hyperscanning_DHMC/sub-pair',...
    num2str(pairN), 'DHMC_fmriprep/fmriprep/sub-pair', num2str(pairN),...
    'DHMC/func/sub-pair', num2str(pairN), 'DHMC_task-storytelling_',...
    'acq-3mm_run-', num2str(runN), '_bold_space-MNI152NLin2009cAsym_', ...
    'preproc_nuisRegr_speaker.mat'];
if ~exist(dhmcFileSpeaker, 'file')
    error([char(10), 'Cannot find dhmcFile at ', dhmcFileSpeaker]);
end
% get parts of filename
[dhmcFileSpeakerPath, dhmcFileSpeakerBase, ~] = fileparts(dhmcFileSpeaker);

dhmcFileListener = ['/flash/wheatley/adamb/hyperscanning_DHMC/sub-pair',...
    num2str(pairN), 'DHMC_fmriprep/fmriprep/sub-pair', num2str(pairN),...
    'DHMC/func/sub-pair', num2str(pairN), 'DHMC_task-storytelling_',...
    'acq-3mm_run-', num2str(runN), '_bold_space-MNI152NLin2009cAsym_', ...
    'preproc_nuisRegr_listener.mat'];
if ~exist(dhmcFileListener, 'file')
    error([char(10), 'Cannot find dhmcFile at ', dhmcFileListener]);
end
% get parts of filename
[dhmcFileListenerPath, dhmcFileListenerBase, ~] = fileparts(dhmcFileListener);

disp([char(10), 'Got all filenames for EPI data mat files']);


%% Save files

dhmcSpeakerSaveF = ['/flash/wheatley/adamb/matlab/DHMC_DBIC/pair',num2str(pairN),'/coupling_run',num2str(runN),'_1.mat'];
dbicSpeakerSaveF = ['/flash/wheatley/adamb/matlab/DHMC_DBIC/pair',num2str(pairN),'/coupling_run',num2str(runN),'_2.mat'];


%% Load data

load(dbicFileSpeaker); % results in variable dbicSpeaker
load(dbicFileListener); % results in variable dbicListener
load(dhmcFileSpeaker); % results in variable dhmcSpeaker
load(dhmcFileListener); % results in variable dhmcListener

disp([char(10), 'Loaded speaker-listener parsed data']);


%% Call couplingFMRI

disp([char(10), 'Calling couplingFMRI on first half of data, dhmcSpeaker and dbicListener']);
if stats
    [b, Rsq, F, pF, pP] = couplingFMRI(dhmcSpeaker, dbicListener, 3, stats);
    save(dhmcSpeakerSaveF, 'b', 'Rsq', 'F', 'pF', 'pP');
    clearvars b Rsq F pF pP;
else
    [b, Rsq, F, pF] = couplingFMRI(dhmcSpeaker, dbicListener, 3);
    save(dhmcSpeakerSaveF, 'b', 'Rsq', 'F', 'pF');
    clearvars b Rsq F pF;
end

disp([char(10), 'Calling couplingFMRI on second half of data, dbicSpeaker and dhmcListener']);
if stats
    [b, Rsq, F, pF, pP] = couplingFMRI(dbicSpeaker, dhmcListener, 3, stats);
    save(dbicSpeakerSaveF, 'b', 'Rsq', 'F', 'pF', 'pP');
    clearvars b Rsq F pF pP;
else
    [b, Rsq, F, pF] = couplingFMRI(dbicSpeaker, dhmcListener, 3);
    save(dbicSpeakerSaveF, 'b', 'Rsq', 'F', 'pF');
    clearvars b Rsq F pF;
end

disp([char(10), 'Done!']);



return








