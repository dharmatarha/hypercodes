function [packetStats, timeData] = timestampsAnalyze_First(pairN, runN)

%% Function for basic analysis of audio packet lags and losses
%
% !!! First dataset !!!
%
% [packetStats, timeData] = timestampsAnalyze_First(pairN, runN)
%
% Inputs:
% pairN - Pair number, integer between 1 - 14
% runN - Run number, either 1 or 2
%
% Output:
% packetStats - Basic audio lag infor in struct, both for dbic and dhmc
% timeData - struct with all audio lag data, both for dbic and dhmc
%
% Logic:
% - We take the two timestamps.csv files from the two ends.
% - Lags are a result of network transmission time + clock drift, we
% estimate both by looking at the differences in lags across sites
% (assuming clock difference is constant for a given run)
% - We get M, SD for network time and a point estimate of clock drift
%



%% Input checks

if nargin ~= 2
    error('Need input args pairN and runN');
end
if ~ismember(pairN,1:14)
    error('Input arg pairN should be int between 1 and 14');
end
if ismember(pairN, 1:13)
    if ~ismember(runN,1:2)
        error('Input arg runN should be 1 or 2');
    end
elseif pairN == 14
    if ~ismember(runN,[1,2,5])
        error('Input arg runN should be 1 or 2 or 5'); % 5 is for repeated run1
    end
end


%% Basics

% get path to log files
if ismember(pairN, 1:13)
    logPath = getLogPath_First(pairN, runN);
elseif pairN == 14
    logPath = getLogPath_14(pairN, runN);
end
    
% result struct
packetStats = struct;
% temp struct for storing timing data
timeData = struct;
% threshold for lag - any lag above this value triggers a warning
threshold = 0.5;

disp([char(10), 'Started timestampsAnalyze_First script with input args ',...
    num2str(pairN), ' (pair number) and ', num2str(runN), ' (run number)']);
disp(['Lags above ', num2str(threshold), ' sec will be reported']);


%% Read in audio packet timestamp data

for site = {'dbic', 'dhmc'}
    filename = [logPath.(site{1}), 'timestamps.csv'];
    timeData.(site{1}) = csvread(filename);
end

disp([char(10), 'Loaded audio timestamp log files']);


%% Get stats

for site = {'dbic', 'dhmc'}
    % ratio of audio packets successfully transmitted
    packetStats.(site{1}).packetRatio = size(timeData.(site{1}), 1) / timeData.(site{1})(end,1);
    % mean lag
    packetStats.(site{1}).meanLag = mean(timeData.(site{1})(:,4));
    % sd lag
    packetStats.(site{1}).sdLag = std(timeData.(site{1})(:,4));
    % feedback about large lags
    if any(timeData.(site{1})(:,4) > threshold)
        [idx, ~] = find(timeData.(site{1})(:,4) > threshold);
        disp([char(10), 'Found audio packets with above-threshold lags ',...
            'in ', num2str(size(idx,1)), ' cases. Indices: ',]);
        disp(idx);
    end
end

% get network time + clock difference estimates
packetStats.networkTime = (packetStats.dbic.meanLag + packetStats.dhmc.meanLag)/2;
packetStats.clockDiff = abs((packetStats.dbic.meanLag - packetStats.dhmc.meanLag)/2);

disp([char(10), 'Got timestamp info.']);
disp(['Network time estimate: ', num2str(packetStats.networkTime), ' sec']);
disp(['Clok difference estimate: ', num2str(packetStats.clockDiff), ' sec']);


return
