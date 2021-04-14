function parseEPI_listeningTask(pairN,debug)

% no interpolation necessary as we don't have auditory stimulus timestamps
% beyond the start and finish
% timingInfo_runs34.m will output the relevant task-related TR indices
% (including the TR just before audio onset and the first TR following
% audio offset)


% Check inputs
if nargin < 2
    debug = 0; % default to loading data from drzeuss
end


% set folder names and load epi files
if debug
    timeFolder = '~/Documents/MATLAB/Wheatley/';
    epiFolder = timeFolder;
else
    timeFolder = '/afs/.dbic.dartmouth.edu/usr/wheatley/jd/timingInfo_runs34/';
    epiFolder = '/afs/.dbic.dartmouth.edu/usr/wheatley/jd/control_tasks/';
end

% load timingInfo data
load([timeFolder 'timingInfo_' num2str(pairN) '_3.mat']);


% get the DBIC subject name for given pair
if debug
    load('~/Documents/MATLAB/Wheatley/isc_test/pairings2.mat');
else
    load('/flash/wheatley/adamb/matlab/pairings2.mat');
end
cbsSub = pairings2.lookup{pairN, 3}; % e.g. 'hid000002'
dbicSub = pairings2.lookup{pairN, 2}; % e.g. 'sid0000522'

% load EPI timeseries
epiData = cell(2,1); % cell 1 = CBS, cell 2 = DBIC
filenames = cell(2,1);
filenames{1} = [epiFolder 'sub-' cbsSub '_ses-pair0' num2str(pairN) '_task-storytelling3_run-03_bold_space-MNI152NLin2009cAsym_preproc_nuisRegr_2021'];
load([filenames{1} '.mat']) % CBS sub
epiData{1} = tseries;
clear tseries;
filenames{2} = [epiFolder 'sub-' dbicSub '_ses-pair0' num2str(pairN) '_task-storytelling3_run-03_bold_space-MNI152NLin2009cAsym_preproc_nuisRegr_2021'];
load([filenames{2} '.mat']) % DBIC sub
epiData{2} = tseries;
clear tseries;

% approximate number of TRs by which to shift things to use the desired 
% hemodynamic delay
hemoDelay = 6; % hemodynamic delay [s]
trDur = 0.727; % TR duration [s]
TRshift = round(hemoDelay / trDur); % number of TRs by which to shift [TRs]

sites = {'cbs','dbic'};
for SITE = 1:2
    
    % adjust timeseries for hemodynamic delay by removing first TRshift TRs
    epiData{SITE}(1:TRshift,:) = [];
    
    % subset according to task-relevant TRs from timingInfo
    if SITE == 1
        disp([char(10) 'Subsetting EPI data from listening task for CBS sub ' cbsSub ' using hemodynamic delay of ' num2str(TRshift) ' TRs (~' num2str(hemoDelay) ' s)'])
    else
        disp([char(10) 'Subsetting EPI data from listening task for DBIC sub ' dbicSub ' using hemodynamic delay of ' num2str(TRshift) ' TRs (~' num2str(hemoDelay) ' s)'])
    end
    tseries = epiData{SITE}(TTLtaskInds.(sites{SITE}),:);
    
    % delete the nuisance regression .mat file -- unfortunately have to do
    % this due to space restrictions -- but if you should ever need to go
    % back to the nuisance regression .mat file, just rerun the nuisance
    % regression python script (nuisRegr_runs_34_2021.py -- but note that
    % this currently doesn't have an option for processing single subjects
    % -- you'll end up running the nusance regression on everything in
    % whatever folder you select -- can hack this by making a temp folder
    % and just sticking the subject you want's data in there.
    delete([filenames{SITE} '.mat'])
    
    % save
    disp([char(10) 'Saving new timeseries to: ' filenames{SITE} '_interp'])
    save([filenames{SITE} '_interp'],'tseries') 
    
end


