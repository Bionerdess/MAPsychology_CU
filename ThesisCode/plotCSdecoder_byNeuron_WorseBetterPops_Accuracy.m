%% By neuron US decoder accuracy plot pops

clear all

dataDir = 'C:\Users\LabAdmin\Dropbox\ANALYSIS\Anna-Lena\NAccVTA_April2021\'; % combined Data location
saveDir = 'C:\Users\LabAdmin\Dropbox\ANALYSIS\Anna-Lena\NAccVTA_April2021\images\';

pops = load([dataDir,'Populations_byNeuron_WorseBetter.mat']);
worseThan = pops.worseThan;
betterThan = pops.betterThan;
accuracy = load([dataDir,'Accuracy_byNeuron_USdecoder']);
accuracyBL = accuracy.accuracyBL;
accuracy = accuracy.accuracy;

pop1Acc = accuracy(worseThan,:);
pop2Acc = accuracy(betterThan,:);
pop1BL = accuracyBL(worseThan,:);
pop2BL = accuracyBL(betterThan,:);

meanpop1 = mean(pop1Acc,1);
meanpop2 = mean(pop2Acc,1);
sempop1 = std(pop1Acc)/sqrt(size(pop1Acc,1));
sempop2 = std(pop2Acc)/sqrt(size(pop2Acc,1));
meanpop1BL = mean(pop1BL,1);
meanpop2BL = mean(pop2BL,1);
sempop1BL = std(pop1BL)/sqrt(size(pop1BL,1));
sempop2BL = std(pop2BL)/sqrt(size(pop2BL,1));

x = 1:size(meanpop1,2);

fh1 = figure;
figure(fh1)
fill([x fliplr(x)], [meanpop1BL-sempop1BL fliplr(meanpop1BL+sempop1BL)],[191 191 191]/255)
hold on
fill([x fliplr(x)], [meanpop2BL-sempop2BL fliplr(meanpop2BL+sempop2BL)],[191 191 191]/255)
hold on
plot(meanpop1BL,'LineWidth',3)
hold on
plot(meanpop2BL,'LineWidth',3)
hold on
fill([x fliplr(x)], [meanpop1-sempop1 fliplr(meanpop1+sempop1)],[191 191 191]/255)
hold on
fill([x fliplr(x)], [meanpop2-sempop2 fliplr(meanpop2+sempop2)],[191 191 191]/255)
hold on
plot(meanpop1,'LineWidth',3)
hold on
plot(meanpop2,'LineWidth',3)
legend('','','WorseThan BL','BetterThan BL','','','WorseThan','BetterThan')
xlabel('Time [sec]')
ylabel('Accuracy')
ylim([0 100])

fh1.Position = [50 50 800 400];
saveas(fh1,[saveDir,'CSdecoder_byNeuronAccuracy_WorseBetterPops3-8_4.svg'],'svg')

fh2 = figure;
figure(fh2)
fill([x fliplr(x)], [meanpop1BL-sempop1BL fliplr(meanpop1BL+sempop1BL)],[191 191 191]/255)
hold on
plot(meanpop1BL,'LineWidth',3)
hold on
fill([x fliplr(x)], [meanpop1-sempop1 fliplr(meanpop1+sempop1)],[191 191 191]/255)
hold on
plot(meanpop1,'LineWidth',3)
hold on
legend('','WorseThan BL','','WorseThan')
xlabel('Time [sec]')
ylabel('Accuracy')
ylim([0 100])

fh2.Position = [50 50 800 400];
saveas(fh2,[saveDir,'CSdecoder_byNeuronAccuracy_WorsePops3-8_4.svg'],'svg')

fh3 = figure;
figure(fh3)
fill([x fliplr(x)], [meanpop2BL-sempop2BL fliplr(meanpop2BL+sempop2BL)],[191 191 191]/255)
hold on
plot(meanpop2BL,'LineWidth',3)
hold on
fill([x fliplr(x)], [meanpop2-sempop2 fliplr(meanpop2+sempop2)],[191 191 191]/255)
hold on
plot(meanpop2,'LineWidth',3)
hold on
legend('','BetterThan BL','','BetterThan')
xlabel('Time [sec]')
ylabel('Accuracy')
ylim([0 100])

fh3.Position = [50 50 800 400];
saveas(fh3,[saveDir,'CSdecoder_byNeuronAccuracy_BetterPops3-8_4.svg'],'svg')

close all