function confoundsPlot(subject)

% plotting confounds for each subject

baseF = '/home/adamb/Documents/hyperscanning_DBIC/';
load([baseF,'subject',subject,'_1.mat']);
run1 = results;
clear results;
load([baseF,'subject',subject,'_2.mat']);
run2 = results;
clear results;
%run3 = load([baseF,'subject',subject,'_3.mat']);
%run4 = load([baseF,'subject',subject,'_4.mat']);

% loaf the first real frames data
load([baseF,'firstFrames.mat']);
% for DBIC data
if subject=='522'
    firstFrames=firstFrames(1,:);
elseif subject=='120'
    firstFrames=firstFrames(2,:);
elseif subject=='199'
    firstFrames=firstFrames(3,:);
elseif subject == '611'
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
plot([1:1238], run1.FramewiseDisplacement(firstFrames(1,1):firstFrames(1,1)+1237));
title(['FrWiseDispl, sub sid',subject,', story run 1']);
xlabel('TRs');
ylabel('mm');
ylim([-1.5, 1.5]);

print (['subject_sid',subject,'_run1_FWD'],'-djpeg');

close all;

% Movement factors - X, Y, Z
plot([1:1238], run1.X(firstFrames(1,1):firstFrames(1,1)+1237), 'b');
hold on;
plot([1:1238], run1.Y(firstFrames(1,1):firstFrames(1,1)+1237), 'r');
plot([1:1238], run1.Z(firstFrames(1,1):firstFrames(1,1)+1237), 'g');
title(['Motion factors X, Y, Z, sub sid',subject,', story run 1']);
xlabel('TRs');
ylabel('mm');
legend('X', 'Y', 'Z');
ylim([-1.5, 1.5]);
hold off;

print (['subject_sid',subject,'_run1_Motion'],'-djpeg');

close all;

% Movement factors - RotX, RotY, RotZ
plot([1:1238], run1.RotX(firstFrames(1,1):firstFrames(1,1)+1237), 'b');
hold on;
plot([1:1238], run1.RotY(firstFrames(1,1):firstFrames(1,1)+1237), 'r');
plot([1:1238], run1.RotZ(firstFrames(1,1):firstFrames(1,1)+1237), 'g');
title(['Motion factors RotX, RotY, RotZ, sub sid',subject,', story run 1']);
xlabel('TRs');
ylabel('rad');
legend('RotX', 'RotY', 'RotZ');
ylim([-0.06, 0.06]);
hold off;

print (['subject_sid',subject,'_run1_MotionRot'],'-djpeg');

close all;


%% Run 2

% storytelling 2
plot([1:1238], run2.FramewiseDisplacement(firstFrames(1,2):firstFrames(1,2)+1237));
title(['FrWiseDispl, sub sid',subject,', story run 2']);
xlabel('TRs');
ylabel('mm');
ylim([-1.5, 1.5]);

print (['subject_sid',subject,'_run2_FWD'],'-djpeg');

close all;

% Movement factors - X, Y, Z
plot([1:1238], run2.X(firstFrames(1,2):firstFrames(1,2)+1237), 'b');
hold on;
plot([1:1238], run2.Y(firstFrames(1,2):firstFrames(1,2)+1237), 'r');
plot([1:1238], run2.Z(firstFrames(1,2):firstFrames(1,2)+1237), 'g');
title(['Motion factors X, Y, Z, sub sid',subject,', story run 2']);
xlabel('TRs');
ylabel('mm');
legend('X', 'Y', 'Z');
ylim([-1.5, 1.5]);
hold off;

print (['subject_sid',subject,'_run2_Motion'],'-djpeg');

close all;

% Movement factors - RotX, RotY, RotZ
plot([1:1238], run2.RotX(firstFrames(1,2):firstFrames(1,2)+1237), 'b');
hold on;
plot([1:1238], run2.RotY(firstFrames(1,2):firstFrames(1,2)+1237), 'r');
plot([1:1238], run2.RotZ(firstFrames(1,2):firstFrames(1,2)+1237), 'g');
title(['Motion factors RotX, RotY, RotZ, sub sid',subject,', story run 2']);
xlabel('TRs');
ylabel('rad');
legend('RotX', 'RotY', 'RotZ');
ylim([-0.06, 0.06]);
hold off;

print (['subject_sid',subject,'_run2_MotionRot'],'-djpeg');

close all;


