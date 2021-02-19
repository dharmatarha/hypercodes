function parseSpeakersFMRI(file, startInMin, endInMin, startMin, endMin)

% parseSpeakersFMRI(file, startInMin, endInMin, startMin, endMin)
%
% Function to cut out the parts of a long fMRI hyperscanning audio when
% someone actually speaks (as the other half is about listening, so
% silence...)
%
% Mandatory inputs:
% file: wav file to work on
% startInMin: Start time (secs) in each minute for the segment to keep 
% endInMin: End time (secs) in each minute for the segment to keep
% startMin: at which minute speech segments start 
% endMin: at which minute speech segments end
%
% Output is the file containing the relevant segments, named
% 'originalFileName_segments.wav'
%
% Example: Audio is 12:30 long, there is a target speech segment from 2:42
% to 3:14, periodically repeating in each minute until 11:14. Then the
% function is called like this:
% parseSpeakersFMRI(file, 42, 14, 2, 11)
%
%


%% Inputs

% check if file exists
if exist(file, 'file') ~= 2
    error('Cannot find file, does it exist?');
end

% check if startInMin and endInMin make sense
if ~ismember(startInMin,0:59) || ~ismember(endInMin,0:59)
    error('Inputs startInMin and endInMin must be integers between 0 and 59');
end

% read in audio immediately
[audio, fs] = audioread(file);
% audio length in secs
audioL = size(audio,1)/fs;

% check if startMin and endMin make sense
if startMin > ceil(audioL/60) || endMin > ceil(audioL/60)
    error('Inputs startMin and endMin cannot be higher then audio length in mins');
end

% don't ask impossible
if startMin > endMin || (startMin == endMin && startInMin >= endInMin)
    error('Something is fishy, check the inputs again...');
end


%% CUT CUT CUT

% length of segments in secs
segmentL = endInMin-startInMin;
if segmentL < 0
    segmentL = segmentL +60;
end

% cut the beginning and the end
audioTrimmed = audio((startMin*60+startInMin)*fs+1 : (endMin*60+endInMin)*fs);

% init result vector
audioSegments = [];
% go through each minute
for i = 1:ceil(size(audioTrimmed,1)/fs/60)
    segment = audioTrimmed( (i-1)*60*fs+1 : (i-1)*60*fs+segmentL*fs );
    audioSegments = [audioSegments;segment];
end


%% Save results

% create file name
[filepath,name,ext] = fileparts(file);
saveFile = [filepath,'/',name,'_segments',ext];

% write out
audiowrite(saveFile, audioSegments, fs);

return



