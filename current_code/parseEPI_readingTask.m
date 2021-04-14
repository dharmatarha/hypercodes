function parseEPI_readingTask(pairN,debug)

% load TTLs
% load stimulus flips
% load epi timeseries
% add whatever hemodynamic delay factor to the stimulus flips to get the
% target times for interpolation

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
load([timeFolder 'timingInfo_' num2str(pairN) '_4.mat']);


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
filenames{1} = [epiFolder 'sub-' cbsSub '_ses-pair0' num2str(pairN) '_task-storytelling4_run-04_bold_space-MNI152NLin2009cAsym_preproc_nuisRegr_2021'];
load(filenames{1}) % CBS sub
epiData{1} = tseries;
clear tseries;
filenames{2} = [epiFolder 'sub-' dbicSub '_ses-pair0' num2str(pairN) '_task-storytelling4_run-04_bold_space-MNI152NLin2009cAsym_preproc_nuisRegr_2021'];
load(filenames{2}) % DBIC sub
epiData{2} = tseries;
clear tseries;

% approximate number of TRs by which to shift things to use the desired 
% hemodynamic delay
hemoDelay = 6; % hemodynamic delay [s]
trDur = 0.727; % TR duration [s]
TRshift = round(hemoDelay / trDur); % number of TRs by which to shift [TRs]

% implement a maximum TR shift based on the lowest number of available TRs 
% following reading stimulus offset (DBIC sub s000416). This limits us to a
% max hemodynamic delay of 5 TRs (about 3.64 seconds). Not ideal, but
% better than nothing.
if TRshift > 5
    TRshift = 5;
    hemoDelay = round(TRshift * trDur * 100) / 100; %get new hemodynamic delay
end 

% set interpolation method for built-in interp1
interpMethod = 'pchip'; % shape preserving piecewise cubic - does not over/undershoot data as it often happens with 'spline'

% for each site...
sites = {'cbs','dbic'};
for SITE = 1:2
    
    % adjust timeseries for hemodynamic delay by removing first TRshift TRs
    epiData{SITE}(1:TRshift,:) = [];
    
    % remove rows off the end of the timeseries so that it's the same length
    % as the ttls vector so that we can use the interp1 function
    rowsToRemove = size(epiData{SITE},1) - length(ttls.(sites{SITE})); % get # of rows to remove
    if rowsToRemove > 0
        epiData{SITE}((end-rowsToRemove+1):end,:) = []; % remove last rowsToRemove rows
    end
    
    % interpolate to get new timeseries
    if SITE == 1
        disp([char(10) 'Running interpolation on reading task EPI data for CBS sub ' cbsSub ' using hemodynamic delay of ' num2str(TRshift) ' TRs (~' num2str(hemoDelay) ' s)'])
    else
        disp([char(10) 'Running interpolation on reading task EPI data for DBIC sub ' dbicSub ' using hemodynamic delay of ' num2str(TRshift) ' TRs (~' num2str(hemoDelay) ' s)'])
    end
    tseries = interp1(ttls.(sites{SITE}), epiData{SITE}, stimFlips.(sites{SITE}), interpMethod);
    
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

 






