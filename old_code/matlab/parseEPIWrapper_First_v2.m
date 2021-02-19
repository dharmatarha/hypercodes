function parseEPIWrapper_First_v2(pairN, runN)

%% Wrapper function around parseEPI_First.m
%
% !!! First dataset, DBIC & DHMC !!!
%
% !!! To be run on DRZEUSS !!!
%
% parseEPIWrapper_First(pairN, runN)
%
% While parseEPI handles the slicing/concatenating of epi timeseries, this
% wrapper does the file handling, output saving, etc.
%
% Inputs:
% pairN:        Pair number, int between 1-14
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

disp([char(10), 'Started parseEPIWrapper_First_v2.m with inputs ',...
    num2str(pairN), ' (pair no) and ', num2str(runN), ' (run no)']);


%% Get timing and epi files for given pairN and runN

% timingInfo output for given pair and run
timingFile = ['/flash/wheatley/adamb/hyperscanning_data_backup/pair',...
    num2str(pairN),'/timingInfo_', num2str(pairN), '_', num2str(runN), '.mat'];
if ~exist(timingFile, 'file')
    error([char(10), 'Cannot find timingInfo file ', timingFile]);
end

% get the DBIC subject name for given pair
load('/flash/wheatley/adamb/matlab/pairings.mat');
dbicSub = pairings.lookup{pairN, 2}; % e.g. 'sid0000522'

% get path to DBIC epi file
dbicFile = ['/flash/wheatley/adamb/hyperscanning_DBIC/sub-', dbicSub,...
    '_fmriprep/fmriprep/sub-', dbicSub, '/func/',...
    'sub-', dbicSub, '_task-storytelling', num2str(runN), '_run-0',... 
    num2str(runN), '_bold_space-MNI152NLin2009cAsym_preproc_nuisRegr.mat'];
if ~exist(dbicFile, 'file')
    error([char(10), 'Cannot find dbicFile at ', dbicFile]);
end
% get parts of filename
[dbicFilePath, dbicFileBase, ~] = fileparts(dbicFile);

% get path to DHMC epi file
dhmcFile = ['/flash/wheatley/adamb/hyperscanning_DHMC/sub-pair',...
    num2str(pairN), 'DHMC_fmriprep/fmriprep/sub-pair', num2str(pairN),...
    'DHMC/func/sub-pair', num2str(pairN), 'DHMC_task-storytelling_',...
    'acq-3mm_run-', num2str(runN), '_bold_space-MNI152NLin2009cAsym_', ...
    'preproc_nuisRegr.mat'];
if ~exist(dhmcFile, 'file')
    error([char(10), 'Cannot find dhmcFile at ', dhmcFile]);
end
% get parts of filename
[dhmcFilePath, dhmcFileBase, ~] = fileparts(dhmcFile);

disp([char(10), 'Got EPI data file path for DBIC and DHMC, calling parseEPI']);


%% Call parseEPI_First.m

% we use default value for turnInterpTRs argument
[dbicSpeaker, dbicListener, dhmcSpeaker, dhmcListener] = ...
parseEPI_First_v2(dbicFile, dhmcFile, timingFile);

disp([char(10), 'Got speaker-listener parsed data from parseEPI']);


%% Save out results

% save dbicSpeaker
savef = [dbicFilePath, '/', dbicFileBase, '_speaker.mat'];
save(savef, 'dbicSpeaker');

% save dbicListener
savef = [dbicFilePath, '/', dbicFileBase, '_listener.mat'];
save(savef, 'dbicListener');

% save dhmcSpeaker
savef = [dhmcFilePath, '/', dhmcFileBase, '_speaker.mat'];
save(savef, 'dhmcSpeaker');

% save dhmcListener
savef = [dhmcFilePath, '/', dhmcFileBase, '_listener.mat'];
save(savef, 'dhmcListener');

disp([char(10), 'Saved out parsed data into separate speaker-listener files',...
    ' closing shop!']);

return





