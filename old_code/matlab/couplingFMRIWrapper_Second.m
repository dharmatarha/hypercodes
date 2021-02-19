function couplingFMRIWrapper_Second(pairN, runN, stats)

%% Wrapper function around couplingFMRI.m
%
% !!! Second dataset, DBIC + CBS !!!
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
if ~ismember(pairN,1:9)
    error('Input arg pairN should be int between 1 and 9');
end
if ~ismember(runN,1:2)
    error('Input arg runN should be 1 or 2');
end

disp([char(10), 'Started couplingFMRIWrapper.m with inputs ',...
    num2str(pairN), ' (pair no) and ', num2str(runN), ' (run no)']);


%% Get file path to everything we need

% get the DBIC subject name for given pair
load('/flash/wheatley/adamb/matlab/pairings2.mat');
dbicSub = pairings2.lookup{pairN, 2}; % e.g. 'sid0000522'

% get path to DBIC epi files
dbicFileSpeaker = ['/flash/wheatley/adamb/hyperscanning_DBIC_ses2/sub-', dbicSub,...
    '_fmriprep/fmriprep/sub-', dbicSub, '/ses-pair0', num2str(pairN), '/func/',...
    'sub-', dbicSub, '_ses-pair0', num2str(pairN), '_task-storytelling', num2str(runN), '_run-0',... 
    num2str(runN), '_bold_space-MNI152NLin2009cAsym_preproc_nuisRegr_speaker.mat'];
if ~exist(dbicFileSpeaker, 'file')
    error([char(10), 'Cannot find dbicFileSpeaker at ', dbicFileSpeaker]);
end
% get parts of filename
[dbicFileSpeakerPath, dbicFileSpeakerBase, ~] = fileparts(dbicFileSpeaker);

dbicFileListener = ['/flash/wheatley/adamb/hyperscanning_DBIC_ses2/sub-', dbicSub,...
    '_fmriprep/fmriprep/sub-', dbicSub, '/ses-pair0', num2str(pairN), '/func/',...
    'sub-', dbicSub, '_ses-pair0', num2str(pairN), '_task-storytelling', num2str(runN), '_run-0',... 
    num2str(runN), '_bold_space-MNI152NLin2009cAsym_preproc_nuisRegr_listener.mat'];
if ~exist(dbicFileListener, 'file')
    error([char(10), 'Cannot find dbicFileSpeaker at ', dbicFileListener]);
end
% get parts of filename
[dbicFileListenerPath, dbicFileListenerBase, ~] = fileparts(dbicFileListener);


% get path to CBS epi files
cbsFileSpeaker = ['/flash/wheatley/adamb/hyperscanning_CBS/sub-hid00000',...
    num2str(pairN), '_fmriprep/fmriprep/sub-hid00000', num2str(pairN),...
    '/ses-pair0', num2str(pairN),'/func/sub-hid00000', num2str(pairN), '_ses-pair0',...
    num2str(pairN),'_task-storytelling', num2str(runN), '_',...
    'run-0', num2str(runN), '_bold_space-MNI152NLin2009cAsym_', ...
    'preproc_nuisRegr_speaker.mat'];
if ~exist(cbsFileSpeaker, 'file')
    error([char(10), 'Cannot find cbsFile at ', cbsFileSpeaker]);
end
% get parts of filename
[cbsFileSpeakerPath, cbsFileSpeakerBase, ~] = fileparts(cbsFileSpeaker);

cbsFileListener = ['/flash/wheatley/adamb/hyperscanning_CBS/sub-hid00000',...
    num2str(pairN), '_fmriprep/fmriprep/sub-hid00000', num2str(pairN),...
    '/ses-pair0', num2str(pairN),'/func/sub-hid00000', num2str(pairN), '_ses-pair0',...
    num2str(pairN),'_task-storytelling', num2str(runN), '_',...
    'run-0', num2str(runN), '_bold_space-MNI152NLin2009cAsym_', ...
    'preproc_nuisRegr_listener.mat'];
if ~exist(cbsFileListener, 'file')
    error([char(10), 'Cannot find cbsFile at ', cbsFileListener]);
end
% get parts of filename
[cbsFileListenerPath, cbsFileListenerBase, ~] = fileparts(cbsFileListener);

disp([char(10), 'Got all filenames for EPI data mat files']);


%% Save files

cbsSpeakerSaveF = ['/flash/wheatley/adamb/matlab/CBS_DBIC/pair',num2str(pairN),'/coupling_run',num2str(runN),'_1.mat'];
dbicSpeakerSaveF = ['/flash/wheatley/adamb/matlab/CBS_DBIC/pair',num2str(pairN),'/coupling_run',num2str(runN),'_2.mat'];


%% Load data

load(dbicFileSpeaker); % results in variable dbicSpeaker
load(dbicFileListener); % results in variable dbicListener
load(cbsFileSpeaker); % results in variable cbsSpeaker
load(cbsFileListener); % results in variable cbsListener

disp([char(10), 'Loaded speaker-listener parsed data']);


%% Call couplingFMRI

disp([char(10), 'Calling couplingFMRI on first half of data, cbsSpeaker and dbicListener']);
if stats
    [b, Rsq, F, pF, pP] = couplingFMRI(cbsSpeaker, dbicListener, 6, stats);
    save(cbsSpeakerSaveF, 'b', 'Rsq', 'F', 'pF', 'pP');
    clearvars b Rsq F pF pP;
else
    [b, Rsq, F, pF] = couplingFMRI(cbsSpeaker, dbicListener, 6);
    save(cbsSpeakerSaveF, 'b', 'Rsq', 'F', 'pF');
    clearvars b Rsq F pF;
end

disp([char(10), 'Calling couplingFMRI on second half of data, dbicSpeaker and cbsListener']);
if stats
    [b, Rsq, F, pF, pP] = couplingFMRI(dbicSpeaker, cbsListener, 6, stats);
    save(dbicSpeakerSaveF, 'b', 'Rsq', 'F', 'pF', 'pP');
    clearvars b Rsq F pF pP;
else
    [b, Rsq, F, pF] = couplingFMRI(dbicSpeaker, cbsListener, 6);
    save(dbicSpeakerSaveF, 'b', 'Rsq', 'F', 'pF');
    clearvars b Rsq F pF;
end

disp([char(10), 'Done!']);



return








