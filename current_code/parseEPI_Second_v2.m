function [dbicSpeaker, dbicListener, cbsSpeaker, cbsListener] = ...
parseEPI_Second_v2(dbicFile, cbsFile, timingInfo, turnInterpTRs, hemoDelay, netAdjust)

%% Function to parse epi data from one pair into speaker-listener parts
%
% !!! Second dataset, DBIC + cbs !!!
%
% [dbicSpeaker, dbicListener, cbsSpeaker, cbsListener] = ...
%  parseEPI_Second(dbicFile, cbsFile, timingInfo, turnInterpTRs)
%
% The function takes two epi timeseries from members of one pair and
% creates two new timeseries. Each of the new timeseries only contains
% portions where one member is speaking and the other is listening. 
%
% Input args:
% dbicFile:      Mat file containing the DBIC member's EPI data in a matrix,
%                in a variable called 'tseries'
% cbsFile:       Mat file containing the cbs member's EPI data in a matrix,
%                in a variable called 'tseries'
% timingInfo:    Mat file containing the output of the timingInfo.m script
%                for current pair and run
% turnInterpTRs: Vector of timings for each speech turn for interpolation,
%                relative to speech turn start, in secs. Default value is
%                [1:tr:turnNo]
% hemoDelay:     Desired hemodynamic delay in seconds. Default = 6
% netAdjust:     Boolean for whether or not to add the estimated network
%                transmission time to listener interpolation timestamps.
%
% Outputs:
% dbicSpeaker:  Matrix containing the re-arranged EPI data concatenating
%               the speech turns where the DBIC member was speaking
% dbicListener: Matrix containing the re-arranged EPI data concatenating
%               the speech turns where the DBIC member was speaking
% cbsSpeaker:   Matrix containing the re-arranged EPI data concatenating
%               the speech turns where the cbs member was speaking
% cbsListener:  Matrix containing the re-arranged EPI data concatenating
%               the speech turns where the cbs member was speaking
%
% This function is to be part of the processing pipeline:
% (1) fmriprep
% (2) confoundsTruncate.py
% (3) nuisRegr.py (diff versions of it)
% (4) timingInfo.m (diff versions of it)
% (5) parseEPI.m (diff versions of it)
% (6) yet to come ISC script

%% Input checks

if nargin < 3
    error('Need input args dbicFile, cbsFIle and timingInfo');
end
if ~exist(dbicFile, 'file')
    error('Input arg dbicFile not found - it should be a mat file');
end
if ~exist(cbsFile, 'file')
    error('Input arg cbsFile not found - it should be a mat file');
end
if ~exist(timingInfo, 'file')
    error('Input arg timingInfo not found - it should be a mat file');
end
if nargin < 4
    turnInterpTRs = []; % initialize turn taking interpolation TR array
end
if nargin < 5  
    hemoDelay = 6; % use default hemodynamic delay of 6 seconds
end
if nargin < 6
   netAdjust = true; % default to adjusting by network transmission time
end

disp([char(10), 'Started parseEPI with dbicFile: ', dbicFile,...
    ';', char(10), 'cbsFile: ', cbsFile, ';', char(10),...
    'timingInfo: ', timingInfo]);


%% Basics: params, data loading

% get # of TRs by which to shift to achieve inputted hemodynamic delay
trDur = 0.727; % TR duration [s]

% dbicFile and cbsFile both contain a "tseries" variable - an epi data
% matrix
load(dbicFile);
dbicData = tseries;
clear tseries;
load(cbsFile);
cbsData = tseries;
clear tseries;
disp([char(10), 'Loaded EPI data matrices']);

% sanity checks about shape of data
if size(dbicData, 1) > size(dbicData, 2) || size(cbsData, 1) > size(cbsData, 2)
    error([char(10), 'EPI data dimensions are strange, more rows than columns']);
end

% timings info is in the following variables: 'packetStats', 'ttls', 
% 'speechTurnStamps', 'retellStartStamp' and 'stimOnsetError' 
load(timingInfo);
disp([char(10), 'Loaded timingInfo mat file']);


%% Edit EPI timeseries to match TTL log

% There are usually more TRs than timestamps due to recording constraints,
% so we need to do some slicing so that no. of ttl timestamps and TRs in 
% data is equal (constraint of the interpolation function). Always the TRs 
% at the end are the ones without timestamps.

if ~isequal(size(ttls.dbic, 1), size(dbicData, 1))
    disp([char(10), 'WARNING! Number of TRs in DBIC EPI data does not equal no. of ',... 
        'DBIC ttl timestamps from timingInfo']);
    disp(['There are ',  num2str(size(ttls.dbic, 1)), ' TTL timestamps, and data for ',...
        num2str(size(dbicData, 1)), ' TRs']);
    % we assume there was more data than TRs, otherwise we are in deep
    % shit...
    dbicData = dbicData(1:size(ttls.dbic, 1), :);
    disp('Resized DBIC data matrix to only include data we have timestamps for');
end
if ~isequal(size(ttls.cbs, 1), size(cbsData, 1))
    disp([char(10), 'WARNING! Number of TRs in cbs EPI data does not equal no. of ',... 
        'cbs ttl timestamps from timingInfo']);
    disp(['There are ',  num2str(size(ttls.cbs, 1)), ' TTL timestamps, and data for ',...
        num2str(size(cbsData, 1)), ' TRs']);
    % we assume there was more data than TRs, otherwise we are in deep
    % shit...
    cbsData = cbsData(1:size(ttls.cbs, 1), :);  
    disp('Resized cbs data matrix to only include data we have timestamps for');
end

%% Set interpolation parameters

% set interpolation method for built-in interp1
interpMethod = 'pchip';  % shape preserving piecewise cubic - does not over/undershoot data as it often happens with 'spline'
% number of speech turns to take into account
turnNo = 30;
% turn length in secs
turnL = 30;
% get no. of voxels
voxelN = size(dbicData, 2);

disp([char(10), 'Loaded pre-set parameters']);

% if not supplied, set timestamps for interpolated data in each speech turn,
% relative to speech turn start
if isempty(turnInterpTRs)
    
    % default is to keep TR sample rate and center time points within
    % speech turn window
    speechTurnDur = 30; % speech turn duration [s]
    trsPerSpeechDur = round(speechTurnDur / trDur); % average TRs per speech turn
    dummy = 0:trDur:speechTurnDur;
    startPoint = (speechTurnDur - dummy(trsPerSpeechDur)) / 2;
    turnInterpTRs = startPoint:trDur:speechTurnDur;

    disp([char(10), 'Input arg turnInterpTRs was not supplied, ',...
        'we use the default: ']);
    disp(turnInterpTRs);
else
    disp([char(10), 'Input arg turnInterpTRs is set to: ']);
    disp(turnInterpTRs);
end
% Check turnInterpTRs
if ~isvector(turnInterpTRs) || any(turnInterpTRs > turnL)
    error([char(10), 'There is a problem with turnInterpTRs argument, it ',...
        'is either not a vector or a member is larger than the turn length']);
end

% add hemodynamic delay to turnInterpTRs
turnInterpTRs = turnInterpTRs + hemoDelay;
disp([char(10), 'Adjusting interpolation window by ' num2str(hemoDelay) ' s to account for hemodynamic delay']);

%% Do same-speaker-speech-turn concatenations

% preallocate matrices
dbicSpeaker = zeros(length(turnInterpTRs)*turnNo/2, voxelN);
dbicListener = zeros(length(turnInterpTRs)*turnNo/2, voxelN);
cbsSpeaker = zeros(length(turnInterpTRs)*turnNo/2, voxelN);
cbsListener = zeros(length(turnInterpTRs)*turnNo/2, voxelN);
disp([char(10), 'Preallocated memory for result matrices']);

% create timestamps for all needed TRs in all speech turns 
turnTimeStamps = zeros(length(turnInterpTRs), turnNo);
counter=0;
for turnStart = speechTurnStamps.cbs(1:turnNo)'
    counter=counter+1;
    turnTimeStamps(:,counter) = (turnInterpTRs + turnStart)';
end

% adjust for network transmission delay if selected, otherwise, proceed to
% interpolation
if netAdjust
    
    % initialize array for site-sepcific turn time stamps
    tts = NaN(size(turnTimeStamps,1),size(turnTimeStamps,2),2);
    
    % for each site...
    for SITE = 1:2 % 1=DBIC, 2=CBS (DBIC was always the first listener)
        
        % get columns in which participant at given site was listener
        listenerCols = SITE:2:(28+SITE);
        
        % add network transmission time to the timestamps in those columns
        tts(:,:,SITE) = turnTimeStamps;
        tts(:,listenerCols,SITE) = tts(:,listenerCols,SITE) + packetStats.networkTime;
        
    end
    
    % interpolate DBIC data with adjusted timestamps
    dbicDataTurns = interp1(ttls.dbic, dbicData, tts(:,:,1), interpMethod);
    % interpolate cbs data with adjusted timestamps
    cbsDataTurns = interp1(ttls.cbs, cbsData, tts(:,:,2), interpMethod);
    
else
   
    % interpolate DBIC data
    dbicDataTurns = interp1(ttls.dbic, dbicData, turnTimeStamps, interpMethod);
    % interpolate cbs data
    cbsDataTurns = interp1(ttls.cbs, cbsData, turnTimeStamps, interpMethod);
    
end

% sanity check result sizes
if ~isequal(size(dbicDataTurns), [length(turnInterpTRs), turnNo, voxelN])
    disp([char(10), 'Stg is wrong, dbicDataTurns has the wrong size: ']);
    disp(size(dbicDataTurns));
    endrk
end
if ~isequal(size(cbsDataTurns), [length(turnInterpTRs), turnNo, voxelN])
    disp([char(10), 'Stg is wrong, cbsDataTurns has the wrong size: ']);
    disp(size(cbsDataTurns));
end

% go through each turn, put the interpolated data into the right result
% matrix
for currentTurn = 1:turnNo
    % cbs was the speaker first for the second dataset, in all cases
    if mod(currentTurn, 2) == 1
        cbsSpeaker((ceil(currentTurn/2)-1)*length(turnInterpTRs)+1:ceil(currentTurn/2)*length(turnInterpTRs), :) = squeeze(cbsDataTurns(:,currentTurn,:));
        dbicListener((ceil(currentTurn/2)-1)*length(turnInterpTRs)+1:ceil(currentTurn/2)*length(turnInterpTRs), :) = squeeze(dbicDataTurns(:,currentTurn,:));
    elseif mod(currentTurn, 2) == 0
        cbsListener((currentTurn/2-1)*length(turnInterpTRs)+1:(currentTurn/2)*length(turnInterpTRs), :) = squeeze(cbsDataTurns(:,currentTurn,:));
        dbicSpeaker((currentTurn/2-1)*length(turnInterpTRs)+1:(currentTurn/2)*length(turnInterpTRs), :) = squeeze(dbicDataTurns(:,currentTurn,:));
    end
end


return