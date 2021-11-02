%% identify period of maximum modulation (excited or inhibited)

clear all

% get data dir
dataDir = 'C:\Users\LabAdmin\Dropbox\ANALYSIS\Anna-Lena\NAccVTA_April2021\'; % combined Data location
saveDir = 'C:\Users\LabAdmin\Dropbox\ANALYSIS\Anna-Lena\NAccVTA_April2021\images\';

FRdata = load([dataDir,'FR_stage3']);
NClass = load([dataDir,'VTA_values_2SD_stage3.mat']);
NClass = NClass.VTA_values;
Dind = find(cellfun(@(NClass) strcmp(NClass,'d'),NClass(:,5))==1);

type = {'A_nr','B_nr'};

ModMax = cell(length(type),1);
ModMin = cell(length(type),1);
for f = 1:length(type)
    ModMax{f} = nan(length(Dind),2);
    ModMin{f} = nan(length(Dind),2);
    for d = 1:length(Dind)
        CS = FRdata.(type{f}).VTA.CS_z{Dind(d),1};
        CS = mean(CS,1);
        
        exc = find(CS==max(CS));
        inh = find(CS==min(CS));
        ModMax{f}(d,1) = max(CS);
        ModMax{f}(d,2) = exc(1);
        ModMin{f}(d,1) = min(CS);
        ModMin{f}(d,2) = inh(1);
    end
end

ModMax{1}(:,3) = 0;
ModMax{2}(:,3) = 1;
ModMin{1}(:,3) = 0;
ModMin{2}(:,3) = 1;

Max = cat(1,ModMax{1},ModMax{2});
Min = cat(1,ModMin{1},ModMin{2});

[bMax,loglMax,Hmax,statMax] = coxphfit(Max(:,2),Max(:,3));
[bMin,loglMin,Hmin,statMin] = coxphfit(Min(:,2),Min(:,3));

disp(['Cox Regression - Max FR Timing: p= ',num2str(statMax.p),'; b= ',num2str(statMax.beta)])
disp(['Cox Regression - Min FR Timing: p= ',num2str(statMin.p),'; b= ',num2str(statMin.beta)])

fh1 = figure;
figure(fh1)
ecdf(ModMax{1}(:,2),'function','survivor')
hold on
ecdf(ModMax{2}(:,2),'function','survivor')
xlim([0 100])
title('Time of Maximum FR')
xlabel('Time')
ylabel('Proportion of Neurons')
legend('A','B')

fh2 = figure;
figure(fh2)
ecdf(ModMin{1}(:,2),'function','survivor')
hold on
ecdf(ModMin{2}(:,2),'function','survivor')
xlim([0 100])
title('Time of Minimum FR')
xlabel('Time')
ylabel('Proportion of Neurons')
legend('A','B')

fh3 = figure;
figure(fh3)
ecdf(ModMax{1}(:,1))
hold on
ecdf(ModMax{2}(:,1))
title('Maximum FR')
xlabel('FR [z-score]')
ylabel('Proportion of Neurons')
legend('A','B')

fh4 = figure;
figure(fh4)
ecdf(ModMin{1}(:,1))
hold on
ecdf(ModMin{2}(:,1))
title('Minimum FR')
xlabel('FR [z-score]')
ylabel('Proportion of Neurons')
legend('A','B')

[hMax,pMax,statsMax] = kstest2(ModMax{1}(:,1),ModMax{2}(:,1));
[hMin,pMin,statsMin] = kstest2(ModMin{1}(:,1),ModMin{2}(:,1));

disp(['KStest - Max FR: p= ',num2str(pMax),'; kstat= ',num2str(statsMax)])
disp(['KStest - Min FR: p= ',num2str(pMin),'; kstat= ',num2str(statsMin)])

saveas(fh1,[saveDir,'MaxModulationTiming_MaxFR_Cox.svg'],'svg')
saveas(fh2,[saveDir,'MaxModulationTiming_MinFR_Cox.svg'],'svg')
saveas(fh3,[saveDir,'MaxModulationFR_MaxFR_CDF.svg'],'svg')
saveas(fh4,[saveDir,'MaxModulationFR_MinFR_CDF.svg'],'svg')

close all