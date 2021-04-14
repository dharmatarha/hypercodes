% estimate_missing_TTL_logs.m
%
% This script "estimates" TTL_logs for DBIC participants who are missing  
% TTL logs for the listening and reading tasks. The resulting TTL logs (one 
% average vector of timestamps for each task), will be saved on drzeuss and
% loaded by the timeline_Second_runs34.m function. The participants to whom
% this pertains are...
%
% Listening task
% S000102 (DBIC, pair 6)
% 
% Reading task
% S000007 (DBIC, pair 2)
% S000009 (DBIC, pair 3)
% S000416 (DBIC, pair 7)

% preallocate TTL and difference cell arrays
targTTLs = cell(1,2);
targTTLs{1} = NaN(495,1); % listening task [TTLs x reference subjects]
targTTLs{2} = NaN(475,3); % reading task [TTLs x reference subjects]
refTTLs = cell(1,2);
refTTLs{1} = NaN(495,7); % listening task [TTLs x reference subjects]
refTTLs{2} = NaN(475,5); % reading task [TTLs x reference subjects]
diffs = cell(1,2);
diffs{1} = NaN(490,7);
diffs{2} = NaN(470,5);
meanDiff = cell(1,2);

% for each task...
for TASK = 3:4 % listening, then reading
   
   % set target and reference pairs
   if TASK == 3 % listening task
      
       targPairs = 6;
       refPairs = [2:5 7:9];
       
   else % reading task
       
       targPairs = [2 3 7];
       refPairs = [4:6 8:9];
       
   end
   
   disp(['Task ' num2str(TASK)])
   
   % for each target pair...
   for PAIR = 1:length(targPairs)
       
       % get paths to DBIC behavioral folders
       logPath = getLogPath_Second_runs34(targPairs(PAIR), TASK);
       
       % get DBIC TTL log
        filename = [logPath.dbic, 'TTL_log.txt']; % standard filename for serial TTL timestamps
        [targTTLs{TASK-2}(1:5,PAIR), ~, ~] = readTTL_log(filename);
   end
   
   % for each reference pair...
   for PAIR = 1:length(refPairs)
       
       % get paths to DBIC behavioral folders
       logPath = getLogPath_Second_runs34(refPairs(PAIR), TASK);
       
       % get DBIC TTL log
        filename = [logPath.dbic, 'TTL_log.txt']; % standard filename for serial TTL timestamps
        [refTTLs{TASK-2}(:,PAIR), ~, ~] = readTTL_log(filename);
   end
   
   % since each target log contains exactly 5 values, we will get the mean
   % difference between TTL 5 and each TTL following it (e.g, TTL 6 - TTL 
   % 5, TTL 7 - TTL 5, ... , TTL 495 - TTL 5) across DBIC subjects. We will 
   % then add these differences to TTL 5 in the target logs to estimate
   % what the full TTL vectors WOULD have been. This is an alternative to
   % the original method of simply adding 0.727 (the TR length) iteratively
   % to TTL 5 to extend out the target TTL logs. The new way should be
   % robust to any systematic shifts in TR timing from the DBIC scanner.
   
   % for each TTL from 6 on...
   for TTL = 6:size(refTTLs{TASK-2},1)
       diffs{TASK-2}(TTL-5,:) = refTTLs{TASK-2}(TTL,:) - refTTLs{TASK-2}(5,:); %get timestamp difference between that TTL and TTL 5
   end
   meanDiff{TASK-2} = mean(diffs{TASK-2},2); %get mean differences across participants
   
   % for each reference pair...
   for PAIR = 1:length(targPairs)
       
       % get estimated TTL timestamps for target pairs
       targTTLs{TASK-2}(6:end,PAIR) = targTTLs{TASK-2}(5,PAIR) + meanDiff{TASK-2};
       
       % save 
       TTLtimes = targTTLs{TASK-2}(:,PAIR);
       save(['/afs/dbic.dartmouth.edu/usr/wheatley/jd/timingInfo_runs34/estimated_DBIC_TTLtimes_run_' num2str(TASK) '_pair_' num2str(targPairs(PAIR))],'TTLtimes')
   end
   
end


