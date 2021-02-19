function listeningLog = readListeningLog(filename)

%% Function to read in a listeningLog.txt file
% for run3 of the hyperscanning TTL logging files,
%
% On DBIC data, this log file contains the important audio start timestamp
% in system time.
%
% On DHMC data, the same file contains all TTL timestamps as well, mixed
% together with all other events.
%
% Input is filepath
%
% Outputs listeningLog is a cell array, first column is the string, second
% is the timstamp, at least for almost all lines of the txt file.
% On DBIC data, cell (7,2) contains the audio start timestamp in system 
% time. (8,2) is the audio end time.
% On DHMC data, things are a bit different because of the TTL intermix, 
% audio start should be at (9,2).
% A more general and robust solution though is to search for "Audio start
% time" and "Audio end time" in the first column of the array. 
%

%% Initialize variables.
delimiter = {',',':'};
% startRow = 6;
% endRow = 6;

%% Format string for each line of text:
%   column1: text (%s)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%s%s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
% textscan(fileID, '%[^\n\r]', startRow-1, 'ReturnOnError', false);
% textscan(fileID, '%[^\n\r]', 'ReturnOnError', false);
% dataArray = textscan(fileID, formatSpec, endRow-startRow+1, 'Delimiter', delimiter, 'ReturnOnError', false);
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Create output variable
listeningLog = [dataArray{1:end}];

return

%% Clear temporary variables
% clearvars filename delimiter startRow endRow formatSpec fileID dataArray ans;