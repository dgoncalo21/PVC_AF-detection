clc, close all, clear all

load('08219_episode-4.mat');
ecg = DAT.ecg;
trueLabeler = DAT.class;
startEnding = DAT.annot;

fs = 250;
windowTime = 10;
windowAF = 250*windowTime;

%% detectar R peaks

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

bpm =(length(peak2)-1)/(length(ecg)/fs)*60

figure
plot (ecg);
hold on
scatter(pos2,peak2)
hold off



%% Checking the signal between R peaks (pos2)

featureVec = [];

for j=3:length(pos2)
   signal = pos2(j-1):pos2(j);
   der = diff(ecg(signal));
   derSum = sum(abs(der)); 
   featureVec = [featureVec, derSum];
end

mediaVec = mean(featureVec);
featureVec(featureVec == 0) = mediaVec*1.6;
mediaVec = mediaVec*1.1;

figure
hold on
stem(featureVec)
plot(1:length(featureVec),ones(length(featureVec))*mediaVec, 'r');
hold off

overThreshold = featureVec > mediaVec;

windowSize = 10;

close all


target = [];
for k=1+windowSize:windowSize/2:length(featureVec)
    scope = overThreshold(k-windowSize:k);
    yesToNoRatio = sum(scope)/length(scope);
    
    if yesToNoRatio > 0.4
        target = [target, ones(1,5)];
    else
        target = [target, zeros(1,5)];
    end
end
target(length(target)+1:length(featureVec)) = 0;


ecgLabel = [];
for j=3:length(pos2)
       signal = pos2(j-1):pos2(j);
       value = target(j-2);
       if value == 0
          ecgLabel = [ecgLabel, zeros(1,length(signal)-1)]; 
       else
          ecgLabel = [ecgLabel, ones(1,length(signal)-1)]; 
       end
end

ecgLabel(length(ecgLabel)+1:length(trueLabeler)) = 0;
ecgLabelFinal = [];

for j=windowAF+1:windowAF/2:length(ecg)
   somaAF = sum(ecgLabel(j-windowAF:j));
   ratio = somaAF/length(j-windowAF:j);
   
   sizer = (length(j-windowAF:j)-1)/2;
   
   if ratio > 0.4
       ecgLabelFinal = [ecgLabelFinal, ones(1,sizer)];
   else
       ecgLabelFinal = [ecgLabelFinal, zeros(1,sizer)];
   end
end

ecgLabelFinal(length(ecgLabelFinal)+1:length(trueLabeler)) = 0;

[confusionMatrix, order] = confusionmat(trueLabeler, ecgLabel);

confusionMatrix

acc = (confusionMatrix(1,1)+confusionMatrix(2,2))/(confusionMatrix(1,1)+confusionMatrix(1,2)+confusionMatrix(2,1)+confusionMatrix(2,2));
sensitivity = confusionMatrix(1,1)/(confusionMatrix(1,1)+confusionMatrix(2,1));
spec = confusionMatrix(2,2)/(confusionMatrix(2,2)+confusionMatrix(1,2));

acc
sensitivity
spec


[confusionMatrix2, order2] = confusionmat(trueLabeler, ecgLabelFinal);

confusionMatrix2

acc2 = (confusionMatrix2(1,1)+confusionMatrix2(2,2))/(confusionMatrix2(1,1)+confusionMatrix2(1,2)+confusionMatrix2(2,1)+confusionMatrix2(2,2));
sensitivity2 = confusionMatrix2(1,1)/(confusionMatrix2(1,1)+confusionMatrix2(2,1));
spec2 = confusionMatrix2(2,2)/(confusionMatrix2(2,2)+confusionMatrix2(1,2));

acc2
sensitivity2
spec2

figure
plot(ecg)
hold on
stem(ecgLabelFinal*0.3)
stem(trueLabeler*-0.3)