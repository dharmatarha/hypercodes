function logPath = getLogPath_Second_runs34(pairN, runN)

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


%% Get path

%initialize
logPath = struct;

%task names
taskNames = {'listening','reading'};

% get DBIC path to logs
logPath.dbic = ['/flash/wheatley/adamb/hyperscanning_data_backup/Harvard_DBIC/',...
    'pair', num2str(pairN), '/DBIC/behav/run', num2str(runN),'_', taskNames{runN-2}, 'Task/'];

% get CBS path to logs
logPath.cbs = ['/flash/wheatley/adamb/hyperscanning_data_backup/Harvard_DBIC/',...
    'pair', num2str(pairN), '/CBS/behav/run', num2str(runN),'_', taskNames{runN-2}, 'Task/'];

    
return