function confoundsPlot_DHMC(subject)

% plotting confounds for each subject

baseF = '/home/adamb/Documents/hyperscanning_DBIC/';
load([baseF,'subject',num2str(subject),'DHMC_1.mat']);
run1 = results; clear results;
load([baseF,'subject',num2str(subject),'DHMC_2.mat']);
run2 = results; clear results;
% run3 = load([baseF,'subject',num2str(subject),'DHMC_3.mat']);
% run4 = load([baseF,'subject',num2str(subject),'DHMC_4.mat']);

% loaf the first real frames data
load([baseF,'firstFrames_DHMC.mat']);
% for DBIC data
if subject==1
    firstFrames=firstFrames(1,:);
elseif subject==2
    firstFrames=firstFrames(2,:);
elseif subject==3
    firstFrames=firstFrames(3,:);
elseif subject == 4
    firstFrames=firstFrames(4,:);
end

% change cell arrays into numbers
a=run1.FramewiseDisplacement;
run1.FramewiseDisplacement = zeros(size(a,1),1);
run1.FramewiseDisplacement(1) = 0;
for i=2:size(a,1)
    run1.FramewiseDisplacement(i) = str2num(a(i,:));
end

a=run2.FramewiseDisplacement;
run2.FramewiseDisplacement = zeros(size(a,1),1);
run2.FramewiseDisplacement(1) = 0;
for i=2:size(a,1)
    run2.FramewiseDisplacement(i) = str2num(a(i,:));
end


%% Run 1

% storytelling 1
plot([1:474], run1.FramewiseDisplacement(firstFrames(1,1):firstFrames(1,1)+473));
title(['FrWiseDispl, sub DHMC',num2str(subject),', story run 1']);
xlabel('TRs');
ylabel('mm');
ylim([-1.5, 1.5]);

print (['subject_DHMC',num2str(subject),'_run1_FWD'],'-djpeg');

close all;

% Movement factors - X, Y, Z
plot([1:474], run1.X(firstFrames(1,1):firstFrames(1,1)+473), 'b');
hold on;
plot([1:474], run1.Y(firstFrames(1,1):firstFrames(1,1)+473), 'r');
plot([1:474], run1.Z(firstFrames(1,1):firstFrames(1,1)+473), 'g');
title(['Motion factors X, Y, Z, sub DHMC',num2str(subject),', story run 1']);
xlabel('TRs');
ylabel('mm');
legend('X', 'Y', 'Z');
ylim([-1.5, 1.5]);
hold off;

print (['subject_DHMC',num2str(subject),'_run1_Motion'],'-djpeg');

close all;

% Movement factors - RotX, RotY, RotZ
plot([1:474], run1.RotX(firstFrames(1,1):firstFrames(1,1)+473), 'b');
hold on;
plot([1:474], run1.RotY(firstFrames(1,1):firstFrames(1,1)+473), 'r');
plot([1:474], run1.RotZ(firstFrames(1,1):firstFrames(1,1)+473), 'g');
title(['Motion factors RotX, RotY, RotZ, sub DHMC',num2str(subject),', story run 1']);
xlabel('TRs');
ylabel('rad');
legend('RotX', 'RotY', 'RotZ');
ylim([-0.06, 0.06]);
hold off;

print (['subject_DHMC',num2str(subject),'_run1_MotionRot'],'-djpeg');

close all;


%% Run 2

% storytelling 2
plot([1:474], run2.FramewiseDisplacement(firstFrames(1,2):firstFrames(1,2)+473));
title(['FrWiseDispl, sub DHMC',num2str(subject),', story run 2']);
xlabel('TRs');
ylabel('mm');
ylim([-1.5, 1.5]);

print (['subject_DHMC',num2str(subject),'_run2_FWD'],'-djpeg');

close all;

% Movement factors - X, Y, Z
plot([1:474], run2.X(firstFrames(1,2):firstFrames(1,2)+473), 'b');
hold on;
plot([1:474], run2.Y(firstFrames(1,2):firstFrames(1,2)+473), 'r');
plot([1:474], run2.Z(firstFrames(1,2):firstFrames(1,2)+473), 'g');
title(['Motion factors X, Y, Z, sub DHMC',num2str(subject),', story run 2']);
xlabel('TRs');
ylabel('mm');
legend('X', 'Y', 'Z');
ylim([-1.5, 1.5]);
hold off;

print (['subject_DHMC',num2str(subject),'_run2_Motion'],'-djpeg');

close all;

% Movement factors - RotX, RotY, RotZ
plot([1:474], run2.RotX(firstFrames(1,2):firstFrames(1,2)+473), 'b');
hold on;
plot([1:474], run2.RotY(firstFrames(1,2):firstFrames(1,2)+473), 'r');
plot([1:474], run2.RotZ(firstFrames(1,2):firstFrames(1,2)+473), 'g');
title(['Motion factors RotX, RotY, RotZ, sub DHMC',num2str(subject),', story run 2']);
xlabel('TRs');
ylabel('rad');
legend('RotX', 'RotY', 'RotZ');
ylim([-0.06, 0.06]);
hold off;

print (['subject_DHMC',num2str(subject),'_run2_MotionRot'],'-djpeg');

close all;


