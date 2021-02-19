function timingInfo_Second(pairN, runN)

%% Wrapper for all timestamp-related functions, gets data for resampling
%
% !!! Second dataset, DBIC + CBS !!!
%
% timingInfo_Second(pairN, runN)
%
% The function call timestampsAnalyze to get packetStats, then call
% timeline to gather all ttl and speech turn timestamps
% Finally, all timing info is rearranged relative to the first event of CBS
% data, meaning that the first speech turn stimulus onset is set to had 
% happened at 0.
% 
% Inputs:
% pairN - pair number, integer between 1 - 14
% runN - run number, only 1 or 2
%
% Outputs are saved out in a struct to the pair level of the behavioral
% data folder tree
%
% Notes:
% - audio buffer length is assumed to be 0.256 sec, but not taken into
% account!!!
% - only run 1-2
% - we correct for clock drift
% 
%

%% Input checks

if nargin ~= 2
    error('Need input args pairN and runN');
end
if ~ismember(pairN,1:9)
    error('Input arg pairN should be int between 1 and 9');
end
if ~ismember(runN,1:2)
    error('Input arg runN should be 1 or 2');
end


%% Basics

% buffer length in secs
bufferL = 0.256;
% path handling
logPath = getLogPath_Second(pairN, runN);
folderParts = regexp(logPath.cbs, 'CBS', 'split');
% savefile
saveFile = [folderParts{1},'timingInfo_', num2str(pairN), '_', num2str(runN), '.mat'];

disp([char(10), 'Started timingInfo_Second script with inputs ', num2str(pairN),...
    ' (pair number) and ', num2str(runN), ' (run number).']);


%% Get network data plus timestamps

% get network transmission estimates
[packetStats, ~] = timestampsAnalyze_Second(pairN, runN);
% get directional clock difference, not absolute value, plus value means
% that DBIC clock was further ahead
if packetStats.dbic.meanLag > packetStats.cbs.meanLag
    packetStats.dbicClock = packetStats.clockDiff;
else
    packetStats.dbicClock = -packetStats.clockDiff;
end

% get event and ttl timestamps
[ttls, speechTurnStamps, retellStartStamp] = timeline_Second(pairN, runN);

disp([char(10), 'Got network and clock estimates, also event and ttl timestamps.']);


%% Adjust timestamps

% get first cbs stim onset
cbsStart = speechTurnStamps.cbs(1);

% adjust all data relative to cbs stim onset
for site = {'dbic', 'cbs'}
    ttls.(site{1}) = ttls.(site{1}) - cbsStart;
    speechTurnStamps.(site{1}) = speechTurnStamps.(site{1}) - cbsStart;
    retellStartStamp.(site{1}) = retellStartStamp.(site{1}) - cbsStart;
end

% correct for average clock difference
ttls.dbic = ttls.dbic - packetStats.dbicClock;
speechTurnStamps.dbic = speechTurnStamps.dbic - packetStats.dbicClock;
retellStartStamp.dbic = retellStartStamp.dbic - packetStats.dbicClock;

disp([char(10), 'Adjusted timing info with cbs start and clock difference']);

% Average stim onset error across sites
stimOnsetError = mean(speechTurnStamps.dbic - speechTurnStamps.cbs);
disp([char(10), 'Average stim onset error across sites: ', num2str(stimOnsetError)]);


%% Save, close

save(saveFile, 'packetStats', 'ttls', 'speechTurnStamps',... 
    'retellStartStamp', 'stimOnsetError');
disp([char(10), 'Saved out results into: ', char(10), saveFile]);


return


