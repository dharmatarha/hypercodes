% parseEPI_Second_v2_wrapper

% set epi directory
epiDir = '/dartfs-hpc/rc/home/z/f00589z/hyperscanning/storytelling/nuisRegr_output_files/';

% set timing file directory
timingDir = '/dartfs-hpc/rc/home/z/f00589z/hyperscanning/storytelling/timingInfo_runs12/';

% set output directory
outputDir = '/dartfs-hpc/rc/home/z/f00589z/hyperscanning/storytelling/parseEPI_output_files/';

% get subject list
subList = readtable('/dartfs-hpc/rc/home/z/f00589z/hyperscanning/misc/hyperscanning_subject_list.csv');

% get number of pairs
numPairs = length(unique(subList.pairNum));

% set hemoynamic delay [s]
hemoDelay = 6; 

% set whether or not to adjust listener time stamps by network transmission time
netAdjust = true;

% set timepoints at which to interpolate on each turn
trDur = 0.727; % TR [s]
speechTurnDur = 30; % speech turn duration [s]
trsPerSpeechDur = round(speechTurnDur / trDur); % average TRs per speech turn
dummy = 0:trDur:speechTurnDur;
startPoint = (speechTurnDur - dummy(trsPerSpeechDur)) / 2;
turnInterpTRs = startPoint:trDur:speechTurnDur;

% save outputs?
saveOutput = false;

% set whether or not to delete nusiance regression output files (may be
% necessary to save space)
deleteNRO = false;

% run parseEPI
for PAIR = 1:numPairs % for each pair

    % get current pair number
    pairNum = subList.pairNum(PAIR);
    
    for RUN = 1:2
        
        % get timing files
        timingInfo = [timingDir 'timingInfo_' num2str(pairNum) '_' num2str(RUN) '.mat'];
    
        % get epi files
        dbicID = ['sub-', subList.subID{PAIR}, '_ses-pair0', num2str(pairNum), '_task-storytelling', num2str(RUN), '_run-0', num2str(RUN), '_bold_space-MNI152NLin2009cAsym_preproc_nuisRegr_2021.mat'];
        dbicFile = [epiDir dbicID];
        [~,dbicFileBase,~] = fileparts(dbicFile);
        cbsID = ['sub-' subList.subID{PAIR + numPairs} '_ses-pair0' num2str(pairNum) '_task-storytelling', num2str(RUN), '_run-0' num2str(RUN) '_bold_space-MNI152NLin2009cAsym_preproc_nuisRegr_2021.mat'];
        cbsFile = [epiDir cbsID];
        [~,cbsFileBase,~] = fileparts(cbsFile);
        
        % run parseEPI
        [dbicSpeaker, dbicListener, cbsSpeaker, cbsListener] = ...
            parseEPI_Second_v3(dbicFile, cbsFile, timingInfo, turnInterpTRs, hemoDelay, netAdjust);
        
        %% Save out results
        if saveOutput

            % save dbicSpeaker
            savef = [outputDir dbicFileBase '_speaker.mat'];
            save(savef, 'dbicSpeaker');

            % save dbicListener
            savef = [outputDir dbicFileBase, '_listener.mat'];
            save(savef, 'dbicListener');

            % save cbsSpeaker
            savef = [outputDir cbsFileBase, '_speaker.mat'];
            save(savef, 'cbsSpeaker');

            % save cbsListener
            savef = [outputDir cbsFileBase, '_listener.mat'];
            save(savef, 'cbsListener');

            disp([char(10), 'Saved out parsed data into separate speaker-listener files',...
                ' closing shop!']);
        end
        
        %% Delete nuisance regression output files if selected
        if deleteNRO
            delete(dbicFile)
            delete(cbsFile)
        end
    end
end