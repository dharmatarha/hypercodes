function [croppedConfs, interpConfs, ttlsTask] = cropConfounds(pairN, runN)

%% Function to crop and resample confounds to task start, for both sites
%
% [croppedConfs, interpConfs, ttlsTask] = cropConfounds(pairN, runN)
%
% Inputs:
% pairN: pair number, int, 1-14
% runN: run number, int, 1-4
%
% Output is a struct of confounds, cropped and resampled so that the first
% data point corresponds to the moment of task start
%
% Relies on all the data-reading helper functions and the established file
% structure of hyperscanning_data_backup
%
% Saves out results into hardcoded folder (see baseFolder in Params
% section)
%
%
% Needs to be extended to runs 3 and 4!!!!
%


%% Input check

if nargin < 2
    error('Need two input args, pairN and runN');
end

if ~ismember(pairN, 1:14)
    error('Input arg pairN should be int between 1 and 14');
end
if ~ismember(runN, 1:4)
    error('Input arg runN should be int between 1 and 14');
end


%% Params

% cropOffset is used when cropping the confound data columns, to ensure
% that the subsequent interpolation step results in reasonable
% interpolation/extrapolation already around the first ("0") time point
%
% cropOffset is time in secs
%
% cropOffset takes different values for the hyperalignment tasks, as there
% are only three TTL's worth of extra data in the beginning of thos tasks
% (that is, only ~2.16 secs for DBIC)
if ismember(runN, 1:2)
    cropOffset = 5;
elseif ismember(runN, 3:4)
    cropOffset = 2;
end
% tNewStart and tNewStep determine the start point and step size for the
% new interpolated / extrapolated timeseries
tNewStart = 0;
tNewStep = 1;
% Treshold for warning, when making a sanity check comparing the two start 
% timestamps across the two sites. Time in secs
offsetTreshold = 0.075;
% base folder to save out the results to
baseFolder = 'resampledConfs_files/';
% stimuli length for run 3 and 4, in secs, ceiling value
stimLengthRun3 = 354;
stimLengthRun4 = 340;


%% Load confounds data

confounds = struct;
confounds.dbic = readConfounds('dbic', pairN, runN);
confounds.dhmc = readConfounds('dhmc', pairN, runN);

disp([char(10),'Loaded confounds tsv files',char(10)]);


%% Get paths for logging files

logPath = getLogPath(pairN, runN);

disp('Got path to log files');


%% Load TTL timestamps
% TTL timestamp files depend on the run and site, 
% DBIC TTLs were saved into files called SerialTTL_log.txt (for runs 1 and
% 2) or TTL_log.txt (for run 3 and 4)
% DHMC timestamps were saved into TTLtimestamps.txt (for runs 1 and 2) or
% into the listeningLog.txt / readingLog.txt files (runs 3 and 4,
% respectively). For the latter type, TTL timestamps are mixed into other
% types of events.

ttls = struct;
taskEndIdx = struct; % no. of ttls recorded

% first for runs 1-2
if ismember(runN, 1:2)
    
    % first DBIC, with serial TTLs
    filename = [logPath.dbic, 'SerialTTL_log.txt']; % standard filename for serial TTL timestamps
    [TTLtimes, ~, ~] = readTTL_log(filename);
    % extract only the absolute times (Unix time) of TTLs
    startTime = TTLtimes(1, 1);
    ttls.dbic = TTLtimes(3:end) + startTime;
    % get number of ttls
    taskEndIdx.dbic = size(ttls.dbic, 1);

    % then DHMC ttls, from TTLtimestamps files
    filename = [logPath.dhmc, 'TTLtimestamps.txt']; % standard filename for key-press type TTL timestamps
    TTLtimestamps = readTTLtimestamps(filename);
    % extract absolute times (Unix time) of timestamps
    ttls.dhmc = cell2mat(TTLtimestamps(2:end, 1));
    % get number of ttls
    taskEndIdx.dhmc = size(ttls.dhmc, 1);
    
% then for runs 3-4
elseif ismember(runN, 3:4)
    
    % for DBIC, the only difference relative to run1 1-2 is in the filename
    filename = [logPath.dbic, 'TTL_log.txt']; % TTL_log.txt
    [TTLtimes, ~, ~] = readTTL_log(filename);
    % extract only the absolute times (Unix time) of TTLs
    startTime = TTLtimes(1, 1);
    ttls.dbic = TTLtimes(3:end) + startTime;
    % get number of ttls
    taskEndIdx.dbic = size(ttls.dbic, 1); 
    
    % for DHMC, we need to read in the log files with different functions,
    % also depending on exact run no.
    if runN == 3
        filename = [logPath.dhmc, 'listeningLog.txt'];
        eventLog = readListeningLog(filename);
    elseif runN == 4
        filename = [logPath.dhmc, 'readingLog.txt'];
        eventLog = readReadingLog(filename);
    end
    % then we need to filter for TTL events, and take the corresponding
    % timestamps - they are already in Unix time
    ttls.dhmc = str2double(eventLog(ismember(eventLog(:, 1),'TTL'), 2));
    % get number of ttls
    taskEndIdx.dhmc = size(ttls.dhmc, 1);
    
end

disp([char(10),'Loaded TTL timestamps, both sites',char(10)]);


%% Load startStamps

startStamp = struct;

% for runs 1-2
if ismember(runN, 1:2)
    
    % for runs 1-2, the event log files are called TimingsLog.txt and have
    % the same structure
    %
    % we loop through sites
    for site = {'dbic', 'dhmc'}
        filename = [logPath.(site{1}), 'TimingsLog.txt'];
        TimingsLog = readTimingsLog(filename);
        startStamp.(site{1}) = str2double(TimingsLog{6, 5});
    end

    % sanity check: the two stamps should be very close to each other
    offset = abs(startStamp.dbic - startStamp.dhmc);
    if offset > offsetTreshold % arbitrary treshold for warning set at params
        warning('startStamp offset seems high!!!');
    end
    disp([char(10), 'startStamp offset between DBIC and DHMC was ', num2str(offset),...
        ' seconds', char(10)]);

% then for run 3
elseif ismember(runN, 3:4)
    
    % loop through sites
    for site = {'dbic', 'dhmc'}
        % for run 3, we need to read in listeningLog.txt and extract 'Audio
        % start time' as stimulus presentation start
        if runN == 3
            filename = [logPath.(site{1}), 'listeningLog.txt'];    
            eventLog = readListeningLog(filename);  
            startTime = eventLog (ismember(eventLog(:, 1), 'Audio start time'), 2);
        % for run 4, we need to read in readingLog.txt and find 'Line 
        % number 0 displayed' as stimulus presentation start
        elseif runN == 4
            filename = [logPath.(site{1}), 'readingLog.txt'];    
            eventLog = readReadingLog(filename);  
            startTime = eventLog (ismember(eventLog(:, 1), 'Line number 0 displayed'), 2);
        end
        % save startTime into startStamp struct
        startStamp.(site{1}) = str2double(startTime);
    end
    
end
    

%% Get task-start TTL indices
% align TTL timestamps relative to task start

taskStartIdx = struct;
ttlsTask = struct;

% loop through sites
for site = {'dbic', 'dhmc'}
    % get timings relative to startStamp
    ttlsTask.(site{1}) = ttls.(site{1}) - startStamp.(site{1});
    % get the index of the first TTL after task start
    [indices, ~] = find(ttlsTask.(site{1}) >= 0);
    taskStartIdx.(site{1}) = indices(1); 
end

disp([char(10),'Got TTL indices relative to task-start',char(10)]);


%% Crop confounds

% init new variable, get list of field names
croppedConfs = confounds;
fields = intersect(fieldnames(croppedConfs.dbic), fieldnames(croppedConfs.dhmc), 'stable');

% iterate over fields (and sites) and crop each confounds columns to values
% during the actual task
for fieldID = 1:length(fields)
    for site = {'dbic', 'dhmc'}
        croppedConfs.(site{1}).(fields{fieldID}) = croppedConfs.(site{1}).(fields{fieldID})(taskStartIdx.(site{1})-cropOffset:taskEndIdx.(site{1}));
    end
end

disp([char(10),'Cropped confounds columns to task start / ending',char(10)]);


%% Interpolate / extrapolate
% for both sites, we want to have values for the same list of timepoints,
% so we interpolate for dbic data and extrapolate for dhmc data

interpConfs = croppedConfs;

% new timepoints series for interpolation - we use a a common series with
% steps of tNewStep sec will the end of the task
tNew = tNewStart:tNewStep:floor(mean([ttlsTask.dhmc(end),ttlsTask.dbic(end)]));

% iterate over fields (and sites) and interpolate / extrapolate columns to values
% during the actual task, using built-in function interpl with 'spine'
% method
for fieldID = 1:length(fields)
    for site = {'dbic', 'dhmc'}
        interpConfs.(site{1}).(fields{fieldID}) = [interp1(ttlsTask.(site{1})(taskStartIdx.(site{1})-cropOffset:taskEndIdx.(site{1})), interpConfs.(site{1}).(fields{fieldID}), tNew, 'spine')]';
    end
end

disp([char(10),'Interpolated cropped confounds columns',char(10)]);


%% Save out everything

% check if baseFolder needs to be created
if ~exist(baseFolder, 'dir')
    mkdir(baseFolder);
end

% save all structs containing the indices for task start, timestamps,
% cropping / resampling details, resampled confounds timeseries
savefile = [baseFolder, 'resampledConfs_pair', num2str(pairN), '_run', num2str(runN), '.mat'];
save(savefile, 'confounds', 'ttls', 'logPath', 'taskStartIdx', 'taskEndIdx', 'startStamp', 'ttlsTask', 'croppedConfs', 'interpConfs', 'tNew');

disp([char(10), 'Saved out results to ', savefile, char(10)]);


return
