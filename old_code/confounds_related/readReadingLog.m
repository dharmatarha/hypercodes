function readingLog = readReadingLog(filename)

%% Function to read in a readingLog.txt file
% for run4 of the hyperscanning TTL logging files,
%
% On DBIC data, this log file contains the important screen display
% timestamps in system time.
%
% On DHMC data, the same file contains all TTL timestamps as well, mixed
% together with all other events.
%
% Input is filepath
%
% Outputs listeningLog is a cell array, first column is the string, second
% is the timstamp, at least for almost all lines of the txt file.
% On DBIC data, cell (7,2) contains the start timestamp for the first line 
% of the text in system time. 
% On DHMC data, things are a bit different because of the TTL intermix, 
% first display timstamp should be at (9,2).
% A more general and robust solution though is to search for "Line number 0
% displayed" in the first column. "Line number 234 displayed" corresponds 
% to the start time of the last line 
%

%% Initialize variables.
delimiter = {',',':'};

%% Format string for each line of text:
%   column1: text (%s)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%s%s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.

dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Create output variable
readingLog = [dataArray{1:end}];

return

%% Clear temporary variables
% clearvars filename delimiter startRow endRow formatSpec fileID dataArray ans;