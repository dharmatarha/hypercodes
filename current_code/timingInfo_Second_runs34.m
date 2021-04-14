function timingInfo_Second_runs34(pairN, runN)

%% Wrapper for all timestamp-related functions, gets data for resampling
%
% !!! Second dataset, DBIC + CBS !!!
%
% timingInfo_Second(pairN, runN)
%
% The function call timestampsAnalyze to get packetStats, then call
% timeline to gather all ttl and speech turn timestamps
% Finally, all timing info is rearranged relative to the first event of CBS
% data, meaning that the reading or listening task onset is set to had 
% happened at 0.
% 
% Inputs:
% pairN - pair number, integer between 1 - 14
% runN - run number, only 3 or 4
%
% Outputs are saved out in a struct to the pair level of the behavioral
% data folder tree
%
% Notes:
% - audio buffer length is assumed to be 0.256 sec, but not taken into
% account!!!
% - only run 3-4
% - we correct for clock drift
% 
% Adapted from timingInfo_Second(), March, 2021 by JDK

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


%% Basics

%buffer length in secs
bufferL = 0.256;

%path handling
logPath = getLogPath_Second_runs34(pairN, runN);
folderParts = regexp(logPath.cbs, 'CBS', 'split');

%savefile
saveFile = ['/afs/.dbic.dartmouth.edu/usr/wheatley/jd/timingInfo_runs34/timingInfo_', num2str(pairN), '_', num2str(runN), '.mat'];

%feedback
disp([char(10), 'Started timingInfo_Second script with inputs ', num2str(pairN),...
    ' (pair number) and ', num2str(runN), ' (run number).']);

%get timestamps and task-specific TRs
[ttls, startStamp, endStamp, TTLtaskInds, stimFlips] = timeline_Second_runs34(pairN, runN);

%% Save, close
save(saveFile, 'ttls', 'startStamp', 'endStamp', 'TTLtaskInds', 'stimFlips');
disp([char(10), 'Saved out results into: ', char(10), saveFile]);

return


