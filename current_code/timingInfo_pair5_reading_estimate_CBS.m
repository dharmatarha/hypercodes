% get DBIC-CBS pair 5 timing info
% requires its own special script because the CBS participant from pair 5 
% is missing a readingLog.txt file

% This script...
% gets DBIC sub from pair 5 (s000535) timing info the normal way

% since CBS sub h00005 (pair 5) does not have a readingLog.txt file, we
% estimate their (zeroed) ttl and stimflip timestamps by taking the mean
% across all other CBS subs

% deal specifically with pair 5 reading task
pairN = 5;
runN = 4;

%savefile
saveFile = ['/afs/.dbic.dartmouth.edu/usr/wheatley/jd/timingInfo_runs34/timingInfo_', num2str(pairN), '_', num2str(runN), '.mat'];

% ttl timestamps structs
ttls = struct;
taskEndIdx = struct; % no. of ttls recorded
startStamp = struct;
endStamp = struct;
stimFlips = struct;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ESTIMATE CBS ttl and stimflip timestamps from the other CBS subs %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%preallocate
ttlData = NaN(473,7); %[ttls x pairs] 
stimData = NaN(236,7); %[stim flips x pairs]

%get zeroed ttl and stimFlip timestamps from other CBS subs
COL = 1;
for PAIR = [2:4 6:9]
    
    %load timingInfo file
    load(['/afs/.dbic.dartmouth.edu/usr/wheatley/jd/timingInfo_runs34/timingInfo_' num2str(PAIR) '_4']);
    
    ttlData(:,COL) = ttls.cbs;
    stimData(:,COL) = stimFlips.cbs;
        
    COL = COL + 1;
    
end

%get means
ttls.cbs = mean(ttlData,2);
stimFlips.cbs = mean(stimData,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% get DBIC timing info the normal way %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%path handling
logPath = getLogPath_Second_runs34(pairN, runN);

%feedback
disp([char(10), 'Started timingInfo_Second script with inputs ', num2str(pairN),...
    ' (pair number) and ', num2str(runN), ' (run number).']);

%get DBIC TTL log
filename = [logPath.dbic, 'TTL_log.txt']; % standard filename for serial TTL timestamps
disp(filename)
[TTLtimes, ~, ~] = readTTL_log(filename);

%if it's one of the pairs with a truncated DBIC TTL_log for the reading 
%task (pairs 2, 3, and 7), extrapolate some fake timestamps...
trDur = 0.727; %set TR duration [s]

%extract only the absolute times (Unix time) of TTLs
startTime = TTLtimes(1, 1);
ttls.dbic = TTLtimes(3:end) + startTime;

%get number of ttls
taskEndIdx.dbic = size(ttls.dbic, 1);

%some feedback
disp([char(10), 'Got TTL timestamps from DBIC data. Found ',...
    num2str(taskEndIdx.dbic), ' TRs'])

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
stimFlips.dbic = []; %initialize flip vector
for ROW = 1:size(TimingsLog,1)
    if (strcmp(TimingsLog{ROW,1},'Line') & strcmp(TimingsLog{ROW,2},'number')) | (strcmp(TimingsLog{ROW,1},'Last') & strcmp(TimingsLog{ROW,2},'line'))
       stimFlips.dbic = [stimFlips.dbic; str2num(TimingsLog{ROW,5})]; 
    end
end

%get indices of start and end stamps and zero TTL, start, and end timestamps to start stamp 
[ttls.dbic, startStamp.dbic, endStamp.dbic, TTLtaskInds.dbic, stimFlips.dbic] = zeroAndGetTasKInds(ttls.dbic,startStamp.dbic,endStamp.dbic, stimFlips.dbic);

%some feedback
disp([char(10), 'Got TTL timestamps and TR indices from DBIC data. Found ', num2str(length(TTLtaskInds.dbic)), ' listening-task-specific TRs from ' num2str(TTLtaskInds.dbic(1)) ' to ' num2str(TTLtaskInds.dbic(end))])

% save, close
save(saveFile, 'ttls', 'startStamp', 'endStamp', 'TTLtaskInds', 'stimFlips');
disp([char(10), 'Saved out results into: ', char(10), saveFile]);