%% Loading Data

clc, clear all, close all
load ('DPVC_106.mat');
ecg=DAT.ecg;

%% Finding R peaks and BPM by filtering

fs = 250;
order = 4;
wc = 20;
fc = wc / (0.5 * fs);
[b, a]=butter(order, fc);
e1 = filter (b, a, ecg);

order = 4;
wc = 5;
fc = wc / (0.5 * fs);
[b,a] = butter(order, fc,'High');
e2 = filter (b, a, e1);

e3 = diff(e2);
e4 = e3.^2;

timeWindow = 0.2;
b = (1/(timeWindow*fs))*ones (1, (timeWindow*fs));
a = 1;
e5 = filter (b, a, e4);

% figure
% plot ([1:length(ecg)]/fs, ecg);
% figure
% plot ([1:length(e1)]/fs, e1);
% figure
% plot ([1:length(e2)]/fs, e2);
% figure
% plot ([1:length(e3)]/fs, e3);
% figure
% plot ([1:length(e4)]/fs, e4);
% figure
% plot ([1:length(e5)]/fs, e5);

threshold = 0.7 * mean (e5);
Pause = 0.3*fs;

peak=[];
i=1;
pos=[];

while i<length(e5)
    
if e5(i)>=threshold
    peak=[peak, e5(i)];
    pos=[pos i];
    i=i+Pause;
else
    i=i+1;
end
end

pos2=[];
peak2=[];

for j=1:length(peak)
    if j==1 
    window=ecg(1:pos(j));
    [M, I]=max(window);
    peak2=[peak2 M];
    pos2=[pos2 I];  
    elseif j==2 
    window=ecg(pos(j-1):pos(j));
    [M, I]=max(window);
    peak2=[peak2 M];
    pos2=[pos2 I]; 
    else
    window=ecg(pos(j)-0.5*fs:pos(j));
    [M, I]=max(window);
    peak2=[peak2 M];
    pos2=[pos2 I+pos(j)-0.5*fs];
    end
end

% figure
% plot ([1:length(ecg)]/fs, ecg);
% xlim([0 2500]);
% hold on
% scatter(pos2/fs,peak2)
% hold off

bpm=(length(peak2)-1)/(length(ecg)/fs)*60


%%

ind=DAT.ind;
class=DAT.pvc;

PPP=[];
NNN=[];
test2=[];
test1=class(1:length(ind));
test1=test1';
areasT=[];

before=7;
after=6;

for i=1:length(ind)
    
    if (ind(i)<before)
        test1 = test1(2:end);
        continue;
    end
    
    if class(i)==0
       window=ecg(ind(i)-before:ind(i)+after);
       area=sum(window-(ecg(ind(i)-before)+ecg(ind(i)+after))/2);
       NNN=[NNN area];
       areasT=[areasT area];
    else
       window2=ecg(ind(i)-before:ind(i)+after);
       area=sum(window2-(ecg(ind(i)-before)+ecg(ind(i)+after))/2);
       PPP=[PPP area];
       areasT=[areasT area];
    end
end

for i=1:length(areasT)
    
    if areasT(i)>= 1500    
        test2=[test2 0];
    else
        test2=[test2 1];
    end
    
end
numN=0;
denN=0;
numP=0;
denP=0;
for i=1:length(test1)
    
    if test1(i)==0
        if test2(i)==0
        numN=numN+1;
        end
        denN=denN+1;
    elseif test1(i)==1
        if test2(i)==1
        numP=numP+1;
        end
        denP=denP+1;
    end
    
end

meanN=mean(NNN);
meanP=mean(PPP);


a=abs(test1-test2);
acerto=1-sum(a)/length(a);
acertoN=numN/denN;
acertoP=numP/denP;

xN=[];
xP=[];
for i=1:length(ind)
    if class(i)==0
        xN=[xN ind(i)];
    end
end
for i=1:length(ind)
    if class(i)==1
        xP=[xP ind(i)];
    end
end

realN=[];
realP=[];

threshold = (meanP+meanN)/2;
totalAreas = [];
for i=2:length(pos2) %ignorado pos2(1)
	window=ecg(pos2(i)-before:pos2(i)+after);
	area=sum(window-(ecg(pos2(i)-before)+ecg(pos2(i)+after))/2);
    if area>= threshold  
        realN=[realN pos2(i)];
    else
        realP=[realP pos2(i)];
    end
    
    totalAreas = [totalAreas, area];
end

figure
plot ([1:length(ecg)]/fs, ecg);
xlim([600 620])
hold on
scatter(xN/fs,ecg(xN),100,'r');% red NORMAL label
scatter(xP/fs,ecg(xP),100,'g'); %green PVC label
scatter(realN/fs,ecg(realN),100,'*r');% red NORMAL
scatter(realP/fs,ecg(realP),100,'*g'); %green PVC
hold off

%% R R Distance

peakDistance = [];

for i=2:size(ind,1)
    peakDistance = [peakDistance, ind(i)-ind(i-1)];
end

label = DAT.pvc;

indNormal = find(label==0);
indPVC = find(label==1);

normalDistance = [];

for i=2:size(indNormal,1)
    normalDistance = [normalDistance, ind(indNormal(i))-ind(indNormal(i)-1)];
end

pvcDistance = [];

for i=2:size(indPVC,1)
    pvcDistance = [pvcDistance, ind(indPVC(i))-ind(indPVC(i)-1)];
end

meanNormal = round(mean(normalDistance));
meanPVC = round(mean(pvcDistance));

errorMargin = 17;

outputPVC = find(peakDistance<meanPVC+errorMargin);
outputNormal = find(peakDistance>=meanPVC+errorMargin);

outputPVC = outputPVC + 1;
outputNormal = outputNormal + 1;

%% Min Dist Clustering

close all

class = class(2:end);
indLess = ind(2:end);
labelLess = label(2:end);

areasT1 = areasT(2:end);

centerNormal = [meanN; meanNormal];
centerPVC = [meanP; meanPVC];
labelClust = [];
allPoints = [];
labelCompare = [];
trueLabelColor = [];

for i=1:length(areasT1)
   area = areasT1(i);
   pDistance = peakDistance (i);
   point = [area; pDistance];
   distN = norm (point - centerNormal);
   distP = norm (point - centerPVC);
   allPoints = [allPoints, point];
   
   
   if distN < distP
       labelClust = [labelClust; 1 0 0];
       labelCompare = [labelCompare; 0];
   else
       labelClust = [labelClust; 0 0 0];
       labelCompare = [labelCompare; 1];

   end
   
   if labelLess(i) == 0
       trueLabelColor = [trueLabelColor; 1 0 0];
   else
       trueLabelColor = [trueLabelColor; 0 0 0];
   end
end
allPoints = allPoints';

%% Plotting centroids, true labels and predicted ones 
figure
hold on
scatter(allPoints(1:end,1), allPoints(1:end,2),[],trueLabelColor);
scatter(allPoints(1:end,1), allPoints(1:end,2),[],labelClust, 'x');
scatter(centerNormal(1),centerNormal(2), 200, '*g');
scatter(centerPVC(1),centerPVC(2), 200, '*b');
title('Clustering with Area and R distance, O are the true labels, X the predicted ones');
hold off

%% Confusion Matrix & Metrics

[confusionMatrix, order] = confusionmat(class, labelCompare);

confusionMatrix

acc = (confusionMatrix(1,1)+confusionMatrix(2,2))/(confusionMatrix(1,1)+confusionMatrix(1,2)+confusionMatrix(2,1)+confusionMatrix(2,2));
sensitivity = confusionMatrix(1,1)/(confusionMatrix(1,1)+confusionMatrix(2,1));
spec = confusionMatrix(2,2)/(confusionMatrix(2,2)+confusionMatrix(1,2));

acc
sensitivity
spec