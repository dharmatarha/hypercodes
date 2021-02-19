function parseEPIWrapper_Second(pairN, runN)

%% Wrapper function around parseEPI_Second.m
%
% !!! Second dataset, DBIC & CBS !!!
%
% !!! To be run on DRZEUSS !!!
%
% parseEPIWrapper_Second(pairN, runN)
%
% While parseEPI handles the slicing/concatenating of epi timeseries, this
% wrapper does the file handling, output saving, etc.
%
% Inputs:
% pairN:        Pair number, int between 1-9
% runN:         Run number, 1 or 2
%
% Outputs are four matrices saved out into the fmriprep preprocessed data
% folder
%
%


%% Input checks

if nargin ~= 2
    error('Need input args pairN and runN');
end
if ~ismember(pairN,1:9)
    error('Input arg pairN should be int between 1 and 14');
end
if ~ismember(runN,1:2)
    error('Input arg runN should be 1 or 2');
end

disp([char(10), 'Started parseEPIWrapper_Second.m with inputs ',...
    num2str(pairN), ' (pair no) and ', num2str(runN), ' (run no)']);


%% Get timing and epi files for given pairN and runN

% timingInfo output for given pair and run
timingFile = ['/flash/wheatley/adamb/hyperscanning_data_backup/Harvard_DBIC/pair',...
    num2str(pairN),'/timingInfo_', num2str(pairN), '_', num2str(runN), '.mat'];
if ~exist(timingFile, 'file')
    error([char(10), 'Cannot find timingInfo file ', timingFile]);
end

% get the DBIC subject name for given pair
load('/flash/wheatley/adamb/matlab/pairings2.mat');
dbicSub = pairings2.lookup{pairN, 2}; % e.g. 'sid0000522'

% get path to DBIC epi file
dbicFile = ['/flash/wheatley/adamb/hyperscanning_DBIC_ses2/sub-', dbicSub,...
    '_fmriprep/fmriprep/sub-', dbicSub, '/ses-pair0', num2str(pairN), '/func/',...
    'sub-', dbicSub, '_ses-pair0', num2str(pairN), '_task-storytelling', num2str(runN), '_run-0',... 
    num2str(runN), '_bold_space-MNI152NLin2009cAsym_preproc_nuisRegr.mat'];
if ~exist(dbicFile, 'file')
    error([char(10), 'Cannot find dbicFile at ', dbicFile]);
end
% get parts of filename
[dbicFilePath, dbicFileBase, ~] = fileparts(dbicFile);

% get path to CBS epi file
cbsFile = ['/flash/wheatley/adamb/hyperscanning_CBS/sub-hid00000',...
    num2str(pairN), '_fmriprep/fmriprep/sub-hid00000', num2str(pairN),...
    '/ses-pair0', num2str(pairN),'/func/sub-hid00000', num2str(pairN), '_ses-pair0',...
    num2str(pairN),'_task-storytelling', num2str(runN), '_',...
    'run-0', num2str(runN), '_bold_space-MNI152NLin2009cAsym_', ...
    'preproc_nuisRegr.mat'];
if ~exist(cbsFile, 'file')
    error([char(10), 'Cannot find cbsFile at ', cbsFile]);
end
% get parts of filename
[cbsFilePath, cbsFileBase, ~] = fileparts(cbsFile);

disp([char(10), 'Got EPI data file path for DBIC and CBS, calling parseEPI']);


%% Call parseEPI_Second.m

% we use default value for turnInterpTRs argument
[dbicSpeaker, dbicListener, cbsSpeaker, cbsListener] = ...
parseEPI_Second(dbicFile, cbsFile, timingFile);

disp([char(10), 'Got speaker-listener parsed data from parseEPI']);


%% Save out results

% save dbicSpeaker
savef = [dbicFilePath, '/', dbicFileBase, '_speaker.mat'];
save(savef, 'dbicSpeaker');

% save dbicListener
savef = [dbicFilePath, '/', dbicFileBase, '_listener.mat'];
save(savef, 'dbicListener');

% save cbsSpeaker
savef = [cbsFilePath, '/', cbsFileBase, '_speaker.mat'];
save(savef, 'cbsSpeaker');

% save cbsListener
savef = [cbsFilePath, '/', cbsFileBase, '_listener.mat'];
save(savef, 'cbsListener');

disp([char(10), 'Saved out parsed data into separate speaker-listener files',...
    ' closing shop!']);

return





