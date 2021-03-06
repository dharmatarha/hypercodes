function [ttls, speechTurnStamps, retellStartStamp] = timeline(pairN, runN)

%% Function to create common timeline for the two fMRI logs
%
% !!! First dataset, DBIC + DHMC !!!
%
% The function loads the log files containing the timestamps of all events,
% both for DBIC and DHMC. One common, "real" timeline is created for all
% these events, taking into account the difference of computer clocks, the
% network transmission time, the audio buffer length.
%
% Inputs:
% pairN - pair number, integer between 1 - 14
% runN - run number, only 1 or 2
% packetStats - packetStats structure as created by the timestampsAnalyze
% function
% bufferL - audio buffer length in seconds
%
% Outputs:
% 

%% Input checks

if nargin ~= 2
    error('Need input args pairN and runN');
end
if ~ismember(pairN,1:14)
    error('Input arg pairN should be int between 1 and 14');
end
if ~ismember(runN,1:2)
    error('Input arg runN should be 1 or 2');
end


%% Basics

pairN = 4;
runN = 1;

% get path to log files
logPath = getLogPath(pairN, runN);
% ttl timestamps structs
ttls = struct;
taskEndIdx = struct; % no. of ttls recorded
% main event timestamp structures
speechTurnStamps = struct;
speechTurnStamps.dbic = zeros(30,1);
speechTurnStamps.dhmc = zeros(30,1);
retellStartStamp = struct;
% text used for timestamp parsing
speechTurnText = {'local', 'time:'};
retellStartText = {'startStampRetell:'};

disp([char(10), 'Started timeline.m script with input arguments ',... 
    num2str(pairN), ' (pair number) and ', num2str(runN),... 
    ' (run number).']);


%% Get TTL timestamps
% this part is from the cropConfounds function

% first DBIC, with serial TTLs
filename = [logPath.dbic, 'SerialTTL_log.txt']; % standard filename for serial TTL timestamps
[TTLtimes, ~, ~] = readTTL_log(filename);
% extract only the absolute times (Unix time) of TTLs
startTime = TTLtimes(1, 1);
ttls.dbic = TTLtimes(3:end) + startTime;
% get number of ttls
taskEndIdx.dbic = size(ttls.dbic, 1);
% some feedback
disp([char(10), 'Got TTL timestamps from DBIC data. Found ',...
    num2str(taskEndIdx.dbic), ' TRs'])

% then DHMC ttls, from TTLtimestamps files
filename = [logPath.dhmc, 'TTLtimestamps.txt']; % standard filename for key-press type TTL timestamps
TTLtimestamps = readTTLtimestamps(filename);
% extract absolute times (Unix time) of timestamps
ttls.dhmc = cell2mat(TTLtimestamps(2:end, 1));
% get number of ttls
taskEndIdx.dhmc = size(ttls.dhmc, 1);
% some feedback
disp([char(10), 'Got TTL timestamps from DHMC data. Found ',...
    num2str(taskEndIdx.dhmc), ' TRs'])


%% Get event timestamps
% here we need to collect speechTurnStamps and retellStartStamp
% again, reusing code from cropConfounds.m

% we loop through sites
for site = {'dbic', 'dhmc'}
    filename = [logPath.(site{1}), 'TimingsLog.txt'];
    TimingsLog = readTimingsLog(filename);
    % get start time
    speechTurnStamps.(site{1})(1) = str2double(TimingsLog{6, 5});
    % go through the whole timingslog file, get each speech turn start +
    % retell start timestamp
    counter = 1;
    for line = 1:size(TimingsLog, 1)
        if isequal(TimingsLog(line, 1:2), speechTurnText)
            counter = counter+1;
            speechTurnStamps.(site{1})(counter) = str2double(TimingsLog(line, 3));
        elseif isequal(TimingsLog(line, 1), retellStartText)
            retellStartStamp.(site{1}) = cell2mat(TimingsLog(line, 2));
            break;
        end 
    end
    % some feedback, sanity check about number of speech turns
    disp([char(10), 'Collected speech turn timestamps for ', site{1}])
    % we expect 31 timestamps, as there is a start one, and then one at
    % the end of each turn 
    if counter ~= 31
        warning(['THERE WERE ', num2str(counter), ' TIMESTAMPS, SOMETHING IS OFF'])
    end
end


return










