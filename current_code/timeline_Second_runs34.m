function [ttls, startStamp, endStamp, TTLtaskInds, stimFlips] = timeline_Second_runs34(pairN, runN)

%% Function to create common timeline for the two fMRI logs
%
% !!! Second dataset, DBIC + CBS !!!
%
% [ttls, speechTurnStamps, retellStartStamp] = timeline_Second_runs34(pairN, runN)
%
% The function loads the log files containing the timestamps of all events,
% both for DBIC and CBS. It finds timestamps for all recorded TTLs, and
% stimulus onsets and offsets. Importantly, it also finds the indices of
% TTLs/TRs that correspond to stimulus presentation.
%
% Inputs:
% pairN - pair number, integer between 1 - 9
% runN - run number, only 3 or 4
%
% Outputs:
% ttls: vector of timestamps for each TR
% startStamp: timestamp for start of task (either audio start or reading start)
% endStamp: timestamp for end of task
% TTLtaskInds: TRs that correspond to the current task
% stimFlips: zeroed timestamps for stimulus flips

%% Input checks

if nargin ~= 2
    error('Need input args pairN and runN');
end
if ~ismember(pairN,1:9)
    error('Input arg pairN should be int between 1 and 9');
end
if ~ismember(runN,3:4)
    error('Input arg runN should be 3 or 4');
end


%% Initialize some stuff

%get path to log files
logPath = getLogPath_Second_runs34(pairN, runN);

%ttl timestamps structs
ttls = struct;
taskEndIdx = struct; % no. of ttls recorded
startStamp = struct;
endStamp = struct;

disp([char(10), 'Started timeline_Second.m script with input arguments ',... 
    num2str(pairN), ' (pair number) and ', num2str(runN),... 
    ' (run number).']);

%% Get DBIC TTL timestamps
% this part is from the cropConfounds function

% if it's one of the pairs with a truncated DBIC TTL_log, load the
% relevant estimated TTLtime variables (estimated using
% 'estimate_missing_TTL_logs.m'). Otherwise get TTL logs the normal way

if (runN == 3 & pairN == 6) | (runN == 4 & (pairN == 2 | pairN == 3 | pairN == 7))  
    
    % load estimated TTL times
    load(['/afs/dbic.dartmouth.edu/usr/wheatley/jd/timingInfo_runs34/estimated_DBIC_TTLtimes_run_' num2str(runN) '_pair_' num2str(pairN)])
    
else
    
    % otherwise get DBIC TTL log the normal way
    filename = [logPath.dbic, 'TTL_log.txt']; % standard filename for serial TTL timestamps
    disp(filename)
    [TTLtimes, ~, ~] = readTTL_log(filename);

end

% extract only the absolute times (Unix time) of TTLs
startTime = TTLtimes(1, 1);
ttls.dbic = TTLtimes(3:end) + startTime;

% get number of ttls
taskEndIdx.dbic = size(ttls.dbic, 1);

% some feedback
disp([char(10), 'Got TTL timestamps from DBIC data. Found ',...
    num2str(taskEndIdx.dbic), ' TRs'])

% initialize TTL task indices structure
TTLtaskInds = struct;

%% Get event timestamps
%the listeningLog.txt and readingLog.txt files are formatted very
%differently, so we handle them separately below. The dbic and cbs .txt
%files are also formatted differently, so they too get handled separately.

% set TR duration [s]
trDur = 0.727; 
% if getting timestamps from the listening task
if runN == 3
    
    %%%%%%%%%%%%
    %%% DBIC %%%
    %%%%%%%%%%%%
    
    %load DBIC listening log
    TimingsLog = readTimingsLog([logPath.dbic, 'listeningLog.txt']); 
    
    %get task start timestamp
    for ROW = 1:size(TimingsLog,1)
        if strcmp(TimingsLog(ROW,1),'Audio') & strcmp(TimingsLog(ROW,2),'start')
            startStamp.dbic = str2double(TimingsLog{ROW, 3}(6:end)); %ignore the 'time:' part at the beginning of the start time cell
            break
        end
    end
    
    %get task end timestamp
    for ROW = 1:size(TimingsLog,1)
        if strcmp(TimingsLog(ROW,1),'Audio') & strcmp(TimingsLog(ROW,2),'end')
            endStamp.dbic = str2double(TimingsLog{ROW, 3}(6:end)); %ignore the 'time:' part at the beginning of the start time cell
            break
        end
    end
    
    %%% Reset start and end timestamps to correspond to TTLs and not %%%
    %%% stimulus onset/offset (DBIC ONLY!) %%%
    
    %reset task start stamp to TTL that comes immediately BEFORE audio onset
    startStamp.dbic = max(ttls.dbic(find(ttls.dbic < startStamp.dbic)));
    
    %reset task end stamp to TTL that comes immediately BEFORE audio offset
    endStamp.dbic = max(ttls.dbic(find(ttls.dbic < endStamp.dbic)));
    
    %get indices of start and end stamps and zero TTL, start, and end timestamps to start stamp 
    [ttls.dbic, startStamp.dbic, endStamp.dbic, TTLtaskInds.dbic, stimFlips.dbic] = zeroAndGetTaskInds(ttls.dbic,startStamp.dbic,endStamp.dbic);
    
    %some feedback
    disp([char(10), 'Got TTL timestamps and TR indices from DBIC data. Found ', num2str(length(TTLtaskInds.dbic)), ' listening-task-specific TRs from ' num2str(TTLtaskInds.dbic(1)) ' to ' num2str(TTLtaskInds.dbic(end))])
    
    %%%%%%%%%%%
    %%% CBS %%%
    %%%%%%%%%%%
    
    %load CBS listening log
    TimingsLog = readTimingsLog([logPath.cbs, 'listeningLog.txt']);
    
    %get task start timestamp
    for ROW = 1:size(TimingsLog,1) %for each row in TimingsLog...
       if strcmp(TimingsLog(ROW,1),'Audio') & strcmp(TimingsLog(ROW,2),'start')
           startStamp.cbs = TimingsLog{ROW - 1, 2}; %using (ROW - 1) here as ROW corresponds to audio onset, but (ROW - 1) corresponds to the last TR that precedes audio onset
           break
       end
    end
    
    %get task end timestamp
    for ROW = 1:size(TimingsLog,1) %for each row in TimingsLog...
        if strcmp(TimingsLog(ROW,1),'Audio') & strcmp(TimingsLog(ROW,2),'end')
            endStamp.cbs = TimingsLog{ROW - 1, 2}; %using ROW - 1 here as ROW corresponds to audio offset, and ROW - 1 corresponds to the last TR that precedes audio offset
           break
        end
    end
    
    %Get CBS TTL timestamps
    ttls.cbs = []; %initialize vector of TTL row indices
    for ROW = 1:size(TimingsLog,1) %for each row in TimingsLog...
        if strcmp(TimingsLog(ROW,1),'TTL,')
            ttls.cbs = [ttls.cbs; TimingsLog{ROW,2}];
        end
    end
    
    %Sanity check - make sure that the difference between each consecutive
    %TTL timestamps is not longer than 1.5x the length of a TR
    for TTL = 2:length(ttls.cbs)
       if ttls.cbs(TTL) - ttls.cbs(TTL-1) > trDur*1.5 
          warning(['The duration between TTLs ' num2str(TTL-1) ' and ' num2str(TTL) ' is greater than 1.5x the TR length (' num2str(trDur*1.5 ) ' s). Investigate...'])
       end
    end
    
    %get indices of start and end stamps and zero TTL, start, and end timestamps to start stamp 
    [ttls.cbs, startStamp.cbs, endStamp.cbs, TTLtaskInds.cbs, stimFlips.cbs] = zeroAndGetTaskInds(ttls.cbs,startStamp.cbs,endStamp.cbs);
    
    %some feedback
    disp([char(10), 'Got TTL timestamps and TR indices from CBS data. Found ', num2str(length(TTLtaskInds.cbs)), ' listening-task-specific TRs from ' num2str(TTLtaskInds.cbs(1)) ' to ' num2str(TTLtaskInds.cbs(end))])
  
%if getting timestamps from the reading task
elseif runN == 4
    
    %%%%%%%%%%%%
    %%% DBIC %%%
    %%%%%%%%%%%%
    
    %load DBIC reading log
    TimingsLog = readTimingsLog([logPath.dbic, 'readingLog.txt']); 
    
    %get timestamp for reading stimulus onset
    for ROW = 1:size(TimingsLog,1)
        if str2double(TimingsLog(ROW,3)) == 0
            startStamp.dbic = str2num(TimingsLog{ROW, 5});
            break
        end
    end
    
    %get timestamp for reading sitmulus offset
    for ROW = 1:size(TimingsLog,1)
        if str2double(TimingsLog(ROW,3)) == 234
            endStamp.dbic = str2num(TimingsLog{ROW + 1, 5});
            break
        end
    end
    
    %%% Reset start and end timestamps to correspond to TTLs and not %%%
    %%% stimulus onset/offset (DBIC ONLY!) %%%
    
    %reset task start stamp to TTL that comes immediately BEFORE audio onset
    startStamp.dbic = max(ttls.dbic(find(ttls.dbic < startStamp.dbic)));
    
    %reset task end stamp to TTL that comes immediately BEFORE audio offset
    endStamp.dbic = max(ttls.dbic(find(ttls.dbic < endStamp.dbic)));
      
    %get vector of stimulus flip timestamps (reading task only)
    stimFlips = struct;
    stimFlips.dbic = []; %initialize flip vector
    for ROW = 1:size(TimingsLog,1)
        if (strcmp(TimingsLog{ROW,1},'Line') & strcmp(TimingsLog{ROW,2},'number')) | (strcmp(TimingsLog{ROW,1},'Last') & strcmp(TimingsLog{ROW,2},'line'))
           stimFlips.dbic = [stimFlips.dbic; str2num(TimingsLog{ROW,5})]; 
        end
    end
    
    %get indices of start and end stamps and zero TTL, start, and end timestamps to start stamp 
    [ttls.dbic, startStamp.dbic, endStamp.dbic, TTLtaskInds.dbic, stimFlips.dbic] = zeroAndGetTaskInds(ttls.dbic,startStamp.dbic,endStamp.dbic, stimFlips.dbic);
    
    %some feedback
    disp([char(10), 'Got TTL timestamps and TR indices from DBIC data. Found ', num2str(length(TTLtaskInds.dbic)), ' listening-task-specific TRs from ' num2str(TTLtaskInds.dbic(1)) ' to ' num2str(TTLtaskInds.dbic(end))])
  
    %%%%%%%%%%%
    %%% CBS %%%
    %%%%%%%%%%%
    
    %load CBS reading log
    TimingsLog = readTimingsLog([logPath.cbs, 'readingLog.txt']);

    %get last TTL timestamp precedeing reading stimulus onset
    for ROW = 1:size(TimingsLog,1)
        if str2double(TimingsLog(ROW,3)) == 0
            startStamp.cbs = TimingsLog{ROW - 2, 2};
            break
        end
    end
    
    %get last TTL timestamp precedeing reading stimulus offset
    for ROW = 1:size(TimingsLog,1)
        if str2double(TimingsLog(ROW,3)) == 234
            endStamp.cbs = TimingsLog{ROW + 3, 2};
            break
        end
    end

    %Get CBS TTL timestamps
    ttls.cbs = []; %initialize vector of TTL row indices
    for ROW = 1:size(TimingsLog,1) %for each row in TimingsLog...
        if strcmp(TimingsLog(ROW,1),'TTL,')
            ttls.cbs = [ttls.cbs; TimingsLog{ROW,2}];
        end
    end
    
    %get vector of stimulus flip timestamps (reading task only)
    stimFlips.cbs = []; %initialize flip vector
    for ROW = 1:size(TimingsLog,1)
        if (strcmp(TimingsLog{ROW,1},'Line') & strcmp(TimingsLog{ROW,2},'number')) | (strcmp(TimingsLog{ROW,1},'Last') & strcmp(TimingsLog{ROW,2},'line'))
           stimFlips.cbs = [stimFlips.cbs; str2num(TimingsLog{ROW,5})]; 
        end
    end
    
    %get indices of start and end stamps and zero TTL, start, and end timestamps to start stamp 
    [ttls.cbs, startStamp.cbs, endStamp.cbs, TTLtaskInds.cbs, stimFlips.cbs] = zeroAndGetTaskInds(ttls.cbs,startStamp.cbs,endStamp.cbs, stimFlips.cbs);

    %Sanity check - make sure that the difference between each consecutive
    %TTL timestamps is not longer than 1.5x the length of a TR
    for TTL = 2:length(ttls.cbs)
       if ttls.cbs(TTL) - ttls.cbs(TTL-1) > trDur*1.5 
          warning(['The duration between TTLs ' num2str(TTL-1) ' and ' num2str(TTL) ' is greater than 1.5x the TR length (' num2str(trDur*1.5 ) ' s). Investigate...'])
       end
    end
    
    %some feedback
    disp([char(10), 'Got TTL timestamps and TR indices from CBS data. Found ', num2str(length(TTLtaskInds.cbs)), ' listening-task-specific TRs from ' num2str(TTLtaskInds.cbs(1)) ' to ' num2str(TTLtaskInds.cbs(end))]);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Reading task timing sanity check %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %get the difference in timestamps for CBS and DBIC stimulus flips
    flipDiffs = stimFlips.cbs - stimFlips.dbic;
    
    %get the difference in timestamps for CBS and DBIC TTLs
    ttlDiffs = ttls.cbs - ttls.dbic;
    
    %%%%%%%%%%%%
    %%% plot %%%
    %%%%%%%%%%%%
    
    %reading stimulus flips
    figure('color','white');
    subplot(2,1,1); hold on
    bar(1:length(flipDiffs),flipDiffs)
    xlabel('reading stimulus #');
    ylabel('CBS - DBIC (s)');
    title('reading stimulus timestamps','FontSize',16);
    yl = ylim; %get y limits
    
    %TTLs
    subplot(2,1,2); hold on
    bar(1:length(ttlDiffs),ttlDiffs)
    xlabel('TR #');
    ylabel('CBS - DBIC (s)');
    title('TR timestamps','FontSize',16);
    ylim(yl); %set y limits equal to that in the stimulus flips sub plot
    
end

end % function end












