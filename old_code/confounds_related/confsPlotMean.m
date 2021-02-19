function confsPlotMean(site)

%% Simple plotting function for the confounds
%
% Include: XYZ, RotXYZ, FD, and aCompCor 1-3
%


%% Check inputs

if ~any(strcmp(site,{'dbic', 'dhmc'}))
    error('Input arg should be dbic or dhmc');
end


%% Magic numbers, params

TRs = 900;
% site = site{1};
basePath = '/home/adamb/Documents/hyperscanning_Confounds/resampledConfs_files/';

if strcmp(site, 'dbic')
    load([basePath, 'confsAllDBIC.mat']);
else
    load([basePath, 'confsAllDHMC.mat']);
end
    
% pairNs = [2:8,10:12,14];
% runNs = 1:2;
% 
% 
% %% Load confounds
% 
% tNewAll = zeros(TRs,11,2);
% XAll = zeros(TRs,11,2);
% YAll = zeros(TRs,11,2);
% ZAll = zeros(TRs,11,2);
% RotAll = zeros(TRs,11,2);
% RotAll = zeros(TRs,11,2);
% RotAll = zeros(TRs,11,2);
% FDAll = zeros(TRs,11,2);
% 
% counter = 1;
% 
% for pairN = pairNs
%     for runN = runNs
%         filename = [basePath, 'resampledConfs_pair', num2str(pairN), '_run', num2str(runN), '.mat'];
%         load(filename);
%         
%         tNewAll(:,counter,runN) = tNew(1:TRs);
%         XAll(:,counter,runN) = interpConfs.(site).X(1:TRs);
%         YAll(:,counter,runN) = interpConfs.(site).Y(1:TRs);
%         ZAll(:,counter,runN) = interpConfs.(site).Z(1:TRs);
%         RotXAll(:,counter,runN) = interpConfs.(site).RotX(1:TRs);
%         RotYAll(:,counter,runN) = interpConfs.(site).RotY(1:TRs);
%         RotZAll(:,counter,runN) = interpConfs.(site).RotZ(1:TRs);
%         FDAll(:,counter,runN) = interpConfs.(site).FramewiseDisplacement(1:TRs);
%     end
%     counter = counter+1;
% end     
        

%% Mean values

% translation
X = squeeze(mean(mean(abs(XAll),3),2));
Y = squeeze(mean(mean(abs(YAll),3),2));
Z = squeeze(mean(mean(abs(ZAll),3),2));

% rotation
RotX = squeeze(mean(mean(abs(RotXAll),3),2));
RotY = squeeze(mean(mean(abs(RotYAll),3),2));
RotZ = squeeze(mean(mean(abs(RotZAll),3),2));

% FD
FD = squeeze(mean(mean(abs(FDAll),3),2));


%% Create a three-subplot plot

h = subplot(1,3,1);
% subplot('Position',[0.1,0.1,0.25,0.5]);
plot(tNewAll(1:TRs,1,1), X(1:TRs), 'b');
hold on;
plot(tNewAll(1:TRs,1,1), Y(1:TRs), 'r');
plot(tNewAll(1:TRs,1,1), Z(1:TRs), 'g');
hold off;
title({'Mean absolute motion factors XYZ,', site, ''});
xlabel('time (sec)');
ylabel('translation (mm)');
legend('X', 'Y', 'Z');
ylim([-0.1, 1]);
xlim([0, TRs]);


h = subplot(1,3,2);
% subplot('Position',[0.4,0.1,0.25,0.5]);
plot(tNewAll(1:TRs,1,1), radtodeg(RotX(1:TRs)), 'b');
hold on;
plot(tNewAll(1:TRs,1,1), radtodeg(RotY(1:TRs)), 'r');
plot(tNewAll(1:TRs,1,1), radtodeg(RotZ(1:TRs)), 'g');
hold off;
title({'Rotation along axes XYZ,', site, ''});
xlabel('time (sec)');
ylabel('rotation (degree)');
legend('RotX', 'RotY', 'RotZ');
ylim([-0.1, 1]);
xlim([0, TRs]);


h = subplot(1,3,3);

plot(tNewAll(1:TRs), FD(1:TRs), 'r');
title({'FramewiseDisplacement,', site, ''});
xlabel('time (sec)');
ylabel('FD (mm)');
ylim([-0.1, 1]);
xlim([0, TRs]);

% h = subplot(2,2,4);
% 
% plot(tNew(1:TRs), interpConfs.(site).aCompCor00(1:TRs), 'b');
% hold on;
% plot(tNew(1:TRs), interpConfs.(site).aCompCor01(1:TRs), 'r');
% plot(tNew(1:TRs), interpConfs.(site).aCompCor02(1:TRs), 'g');
% hold off;
% title({'Anatomically defined noise components (aCompCor)', [site, ' pair',num2str(pairN),', run ', num2str(runN)], ''});
% xlabel('time (sec)');
% ylabel('component weights');
% legend('Comp1', 'Comp2', 'Comp3');
% ylim([-0.1, 0.1]);
% xlim([0, TRs]);


% %% Save
% set(gcf, 'PaperUnits', 'centimeters');
% set(gcf, 'PaperPosition', [0 0 40 10]);
% savefile = [basePath,'meanConfoundsFig_',site,'.tif'];
% saveas(gcf,savefile)


% %% close all cleanup
% close all;

return


