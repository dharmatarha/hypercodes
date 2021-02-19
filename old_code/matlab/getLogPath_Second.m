function logPath = getLogPath_Second(pairN, runN)

%% Helper function to get the exact path to the log files of a given run
%
% !!! Drzeuss version !!!
%
% !!! Second dataset - DBIC and CBS !!!
%
% logPath = getLogPath_Second(pairN, runN)
%
% Inputs are pair number and run number
%
% Output is a struct, with the .dbic and .dhmc fields containing
% corresponding paths
%

%% Check input

if nargin ~= 2
    error('Need input args pairN and runN');
end
if ~ismember(pairN,1:9)
    error('Input arg pairN should be int between 1 and 14');
end
if ~ismember(runN,1:4)
    error('Input arg runN should be int between 1 and 4');
end


%% Init

logPath = struct;
runName = cell(2,1);


%% Get path

% get DBIC path to logs
runName{1} = dir(['/flash/wheatley/adamb/hyperscanning_data_backup/Harvard_DBIC/',...
    'pair', num2str(pairN), '/DBIC/behav/run', num2str(runN),'*']);
logPath.dbic = ['/flash/wheatley/adamb/hyperscanning_data_backup/Harvard_DBIC/',...
    'pair', num2str(pairN), '/DBIC/behav/', runName{1}.name, '/'];

% get DHMC path to logs
runName{2} = dir(['/flash/wheatley/adamb/hyperscanning_data_backup/Harvard_DBIC/',...
    'pair', num2str(pairN), '/CBS/behav/run', num2str(runN),'*']);
logPath.cbs = ['/flash/wheatley/adamb/hyperscanning_data_backup/Harvard_DBIC/',...
    'pair', num2str(pairN), '/CBS/behav/', runName{2}.name, '/'];

    
return