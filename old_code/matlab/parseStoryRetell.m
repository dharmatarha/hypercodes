function parseStoryRetell(file, cutPoint)

% parseStoryRetell(file, cutPoint)
%
% Function to cut a wav into two parts at the specified time point. 
%
% Mandatory inputs:
% file: path to a wav file
% cutPoint: point in time to parse at, in secs from the beginning
%
% Outputs are two files, named [baseFileName]_part[1 OR 2].wav
%

%% Inputs

% check if file exists
if exist(file, 'file') ~= 2
    error('Cannot find file, does it exist?');
end

% read in audio immediately
[audio, fs] = audioread(file);
% audio length in secs
audioL = size(audio,1)/fs;

% check if cutPoint is okay
if cutPoint >= audioL
    error('Arg cutPoint is larger or equal to the full length of audio');
end


%% DO

% audio segments
part1 = audio(1 : cutPoint*fs);
part2 = audio(cutPoint*fs+1 : size(audio,1));

% filenames
[filepath,name,ext] = fileparts(file);
saveFile1 = [filepath,'/',name,'_part1',ext];
saveFile2 = [filepath,'/',name,'_part2',ext];

% save results
audiowrite(saveFile1, part1, fs);
audiowrite(saveFile2, part2, fs);

return



