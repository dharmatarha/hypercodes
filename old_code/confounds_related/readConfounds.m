function confounds = readConfounds(site, pairN, runN, basePath)

% Function to read in a confounds tsv file from the fmriprep output.
%
% Confounds files are expected to be found in a confounds_files folder,
% change the base path if data is elsewhere.
%
% Also expects a pairings.mat for pair number / data_id lookup to be
% present, with columns pairNumber, DBIC_data_ID and DHMC_data_ID in the
% pairings.lookup cell matrix
%
% Inputs:
% site: string, dbic / dhmc
% pairN: int, between 1-14
% runN: int, between 1-4
%
% Returns a struct with all the confounds variables.


%% Check inputs

if nargin == 3
    basePath = '/home/adamb/Documents/hyperscanning_Confounds/confounds_files/';
end

if nargin < 3
    error('Need three input args: site, pairN and runN');
end

if ~any(strcmp(site,{'dbic', 'dhmc'}))
    error('Input arg should be dbic or dhmc');
end
if ~ismember(pairN, 1:14)
    error('Input arg pairN should be between 1-14');
end
if ~ismember(runN, 1:4)
    error('Input arg runN should be between 1-4');
end
if ~exist(basePath, 'dir')
    error('Optional input arg basePath cannot be found');
end


%% Lookup (site, pairN) to subject ID

load pairings.mat;
dataID = pairings.lookup{pairN, (find(strcmp(site, {'dbic', 'dhmc'})))+1};


%% Get exact file

% get the correct filename based on dataID, site, runN - unfortunately,
% naming conventions are different across sites

% if data is from DBIC
if isequal(site, 'dbic')
    filePath = [basePath, 'sub-', dataID, '_task-storytelling', num2str(runN),...
        '_run-0', num2str(runN), '_bold_confounds.tsv'];
% if data is from DHMC
elseif isequal(site, 'dhmc')
    % task name is different depending on runN
    if ismember(runN, 1:2)
        task = 'storytelling';
    elseif runN == 3
        task = 'listening';
    elseif runN == 4
        task = 'reading';
    end
    filePath = [basePath, 'sub-', dataID, '_task-', task, '_acq-3mm_run-',...
        num2str(runN),'_bold_confounds.tsv'];
end


%% Load confounds

% user feedback
disp(['Reading file ',filePath]);

% read
confounds = tdfread(filePath, '\t');


%% Correct char values

fields = {'stdDVARS', 'non0x2DstdDVARS', 'vx0x2DwisestdDVARS', 'FramewiseDisplacement'};
confN = size(confounds.X,1);

for fieldID = fields
    tmp = zeros(confN,1);
    for i = 1:confN
        tmp(i,1) = str2double(confounds.(fieldID{1})(i,:));
    end
    confounds.(fieldID{1}) = tmp;
end


% user feedback
disp('done');

return
