function logPath = getLogPath(pairN, runN)

%% Helper function to get the exact path to the log files of a given run
%
% logPath = getLogPath(pairN, runN)
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
if ~ismember(pairN,1:14)
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
runName{1} = dir(['/media/adamb/TOSHIBA_EXT/hyperscanning_data_backup/',...
    'pair', num2str(pairN), '/DBIC/behav/run', num2str(runN),'*']);
logPath.dbic = ['/media/adamb/TOSHIBA_EXT/hyperscanning_data_backup/',...
    'pair', num2str(pairN), '/DBIC/behav/', runName{1}.name, '/'];

% get DHMC path to logs
runName{2} = dir(['/media/adamb/TOSHIBA_EXT/hyperscanning_data_backup/',...
    'pair', num2str(pairN), '/DHMC/behav/run', num2str(runN),'*']);
logPath.dhmc = ['/media/adamb/TOSHIBA_EXT/hyperscanning_data_backup/',...
    'pair', num2str(pairN), '/DHMC/behav/', runName{2}.name, '/'];

return