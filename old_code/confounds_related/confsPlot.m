function confsPlot(site, pairN, runN, basePath)

%% Simple plotting function for the confounds
%
% confsPlot(site, pairN, runN, basePath)
%
% Include the following confounds: XYZ, RotXYZ, FD, and aCompCor 1-3
%
% Inputs:
% site = 'dbic' or 'dhmc'
% pairN = int between 1 and 14
% runN = int between 1 and 4
% basePath = optional, default should be fine, but str to folder containing
% resampled confounds files
%
% The resulting plot is saved out into a file named:
% [basePath,'confoundsFig_',site,'_pair',num2str(pairN),'_run',num2str(runN),'.tif']
% where basePath is either given as input or the default
% '/home/adamb/Documents/hyperscanning_Confounds/resampledConfs_files/'
%


%% Check inputs

if nargin == 3
    basePath = '/home/adamb/Documents/hyperscanning_Confounds/resampledConfs_files/';
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


%% Magic numbers, params

% the number of resampled TRs to show depends on the run
if ismember(runN, 1:2)
    TRs = 900; % for storytelling runs we only plot the story but no the retell parts
elseif runN == 3
    TRs = 354; % for listening task, the audio stimuli is ~354 secs long
elseif runN == 4
    TRs = 340; % for the reading task the stimuli were ~340 secs long
end
site = site{1};

%% Load confounds

filename = [basePath, 'resampledConfs_pair', num2str(pairN), '_run', num2str(runN), '.mat'];
load(filename);


%% Create a four-subplot plot

h = subplot(2,2,1);

plot(tNew(1:TRs), interpConfs.(site).X(1:TRs), 'b');
hold on;
plot(tNew(1:TRs), interpConfs.(site).Y(1:TRs), 'r');
plot(tNew(1:TRs), interpConfs.(site).Z(1:TRs), 'g');
hold off;
title({'Motion factors XYZ,', [site, ' pair',num2str(pairN),', run ', num2str(runN)], ''});
xlabel('time (sec)');
ylabel('translation (mm)');
legend('X', 'Y', 'Z');
ylim([-1.5, 1.5]);
xlim([0, TRs]);

h = subplot(2,2,2);

plot(tNew(1:TRs), radtodeg(interpConfs.(site).RotX(1:TRs)), 'b');
hold on;
plot(tNew(1:TRs), radtodeg(interpConfs.(site).RotY(1:TRs)), 'r');
plot(tNew(1:TRs), radtodeg(interpConfs.(site).RotZ(1:TRs)), 'g');
hold off;
title({'Rotation along axes XYZ,', [site, ' pair',num2str(pairN),', run ', num2str(runN)], ''});
xlabel('time (sec)');
ylabel('rotation (degree)');
legend('RotX', 'RotY', 'RotZ');
ylim([-1.5, 1.5]);
xlim([0, TRs]);

h = subplot(2,2,3);

plot(tNew(1:TRs), interpConfs.(site).FramewiseDisplacement(1:TRs), 'r');
title({'FramewiseDisplacement,', [site, ' pair',num2str(pairN),', run ', num2str(runN)], ''});
xlabel('time (sec)');
ylabel('FD (mm)');
ylim([-0.4, 1.2]);
xlim([0, TRs]);

h = subplot(2,2,4);

plot(tNew(1:TRs), interpConfs.(site).aCompCor00(1:TRs), 'b');
hold on;
plot(tNew(1:TRs), interpConfs.(site).aCompCor01(1:TRs), 'r');
plot(tNew(1:TRs), interpConfs.(site).aCompCor02(1:TRs), 'g');
hold off;
title({'Anatomically defined noise components (aCompCor)', [site, ' pair',num2str(pairN),', run ', num2str(runN)], ''});
xlabel('time (sec)');
ylabel('component weights');
legend('Comp1', 'Comp2', 'Comp3');
ylim([-0.1, 0.1]);
xlim([0, TRs]);


%% Save
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 30 20]);
savefile = [basePath,'confoundsFig_',site,'_pair',num2str(pairN),'_run',num2str(runN),'.tif'];
saveas(gcf,savefile)


%% close all cleanup
close all;

return


