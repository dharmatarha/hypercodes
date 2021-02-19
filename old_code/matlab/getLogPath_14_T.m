function logPath = getLogPath_14_T(pairN, runN)

%% Helper function to get the exact path to the log files of a given run
%
% !!! Drzeuss version !!!
%
% !!! First dataset - DBIC and DHMC !!!
%
% !!! Only for pair14 !!!
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
if ~isequal(pairN,14)
    error('Input arg pairN should be 14, this version of the script only works on that pair');
end
if ~ismember(runN,1:5)
    error('Input arg runN should be int between 1 and 5');
end


%% Init

logPath = struct;
runName = cell(2,1);


%% Get path

if ismember(runN, 1:4)
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
% if we want to look into the repeated run1 data
elseif runN == 5
    % get DBIC path to logs
    runName{1} = dir(['/media/adamb/TOSHIBA_EXT/hyperscanning_data_backup/',...
        'pair', num2str(pairN), '/DBIC/behav/repeated_run1*']);
    logPath.dbic = ['/media/adamb/TOSHIBA_EXT/hyperscanning_data_backup/',...
        'pair', num2str(pairN), '/DBIC/behav/', runName{1}.name, '/'];
    % get DHMC path to logs
    runName{2} = dir(['/media/adamb/TOSHIBA_EXT/hyperscanning_data_backup/',...
        'pair', num2str(pairN), '/DHMC/behav/repeated_run1*']);
    logPath.dhmc = ['/media/adamb/TOSHIBA_EXT/hyperscanning_data_backup/',...
        'pair', num2str(pairN), '/DHMC/behav/', runName{2}.name, '/'];    
end


return