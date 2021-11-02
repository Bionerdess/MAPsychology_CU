%% calc diff in z-scores A vs B

clear all

% get data dir
dataDir = 'C:\Users\LabAdmin\Dropbox\ANALYSIS\Anna-Lena\NAccVTA_April2021\'; % combined Data location
saveDir = 'C:\Users\LabAdmin\Dropbox\ANALYSIS\Anna-Lena\NAccVTA_April2021\images\';

FR = load([dataDir,'FR_stage3.mat']);
NClass = load([dataDir,'VTA_values_2SD_stage3']); 
NClass = NClass.VTA_values;
Dind = find(cellfun(@(NClass) strcmp(NClass,'d'),NClass(:,5))==1);

diff = nan(length(Dind),150);
diffUS = nan(length(Dind),45);
A = nan(length(Dind),100);
B = nan(length(Dind),100);

for d = 1:length(Dind)
    CSA = FR.A_nr.VTA.CS_z{Dind(d),1};
    CSB = FR.B_nr.VTA.CS_z{Dind(d),1};
    baseA = FR.A_nr.VTA.baseline_z{Dind(d),1};
    baseB = FR.B_nr.VTA.baseline_z{Dind(d),1};
    USA = FR.A_nr.VTA.US_z{Dind(d),1};
    USB = FR.B_nr.VTA.US_z{Dind(d),1};
    
    dataA = [baseA(:,51:100),CSA];
    dataB = [baseB(:,51:100),CSB];
    diff(d,:) = mean(abs(dataA-dataB),1);
    diffUS(d,:) = mean(abs(USA(:,6:50)-USB(:,1:45)),1);
    
end

meanDiff = mean(diff,1);
semDiff = std(diff,1)/sqrt(size(diff,1));
x = 1:150;

fh1 = figure;
figure(fh1)
fill([x fliplr(x)], [meanDiff-semDiff fliplr(meanDiff+semDiff)],[161 161 161]/255)
hold on
plot(meanDiff,'k','LineWidth',3)
hold on
plot([50 50],[0 8],'k')
hold on

ylabel('delta z-score')
xlabel('Time [100ms bins]')
saveas(fh1,[saveDir,'deltaZScore_PSTH_BLCS.svg'],'svg')

meanUS = mean(diffUS,1);
semUS = std(diffUS,1)/sqrt(size(diffUS,1));
xUS = 1:45;

fh2 = figure;
figure(fh2)
fill([xUS fliplr(xUS)], [meanUS-semUS fliplr(meanUS+semUS)],[161 161 161]/255)
hold on
plot(meanUS,'k','LineWidth',3)
ylim([0 8])
ylabel('delta z-score')
xlabel('Time [100ms bins]')
saveas(fh2,[saveDir,'deltaZScore_PSTH_US.svg'],'svg')