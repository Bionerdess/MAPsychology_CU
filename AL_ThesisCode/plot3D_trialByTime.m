%% 3D heatplot - population avgFR during CS and US periods across trials
% x=time bins
% y=trials
% z=avgFR

clear all

% get data dir
dataDir = 'C:\Users\LabAdmin\Dropbox\ANALYSIS\Anna-Lena\NAccVTA_April2021\'; % combined Data location
saveDir = 'C:\Users\LabAdmin\Dropbox\ANALYSIS\Anna-Lena\NAccVTA_April2021\images\';

FRdata = load([dataDir,'FR_','stage3','.mat']);

NClass = load([dataDir,'VTA_values_2SD_','stage3','.mat']);
NClass = NClass.VTA_values;
Dind = find(cellfun(@(NClass) strcmp(NClass,'d'),NClass(:,5))==1);
type = {'A_nr','B_nr'};

nTrials=8;
data = cell(1,length(type));

fh1 = figure;
fh2=figure;
for f=1:length(type)
    data{f} = nan(nTrials,200);
    for t=1:nTrials
        baseline = nan(length(Dind),100);
        CS = nan(length(Dind),100);
        US = nan(length(Dind),50);
        for d=1:length(Dind)
            
            baseline(d,:) = FRdata.(type{f}).VTA.baseline_z{Dind(d),1}(t,:);
            CS(d,:) = FRdata.(type{f}).VTA.CS_z{Dind(d),1}(t,:);
            US(d,:) = FRdata.(type{f}).VTA.US_z{Dind(d),1}(t,:);
        end
        baseline=rmmissing(baseline);
        CS = rmmissing(CS);
        US = rmmissing(US);
        total = [baseline(:,51:100),CS,US];
        data{f}(t,:) = mean(total,1);
    end
end

y=1:size(data{1},1);
x=-5:0.1:14.9;

figure(fh1)
surf(x,y,data{1})
axis([x(1) x(length(x)) y(1) y(length(y)) -2 7])
view(-20,60)
colormap(jet)
colorbar
ylabel('Trial')
xlabel('Time from CS on [sec]')
title('avgFR population across trials - A')

figure(fh2)
surf(x,y,data{2})
axis([x(1) x(length(x)) y(1) y(length(y)) -2 7])
view(-20,60)
colormap(jet)
colorbar
ylabel('Trial')
xlabel('Time from CS on [sec]')
title('avgFR population across trials - B')

saveas(fh1,[saveDir,'Heatplot_avgFRA_2SD_acrossTrials-20-60.svg'],'svg')
saveas(fh2,[saveDir,'Heatplot_avgFRB_2SD_acrossTrials-20-60.svg'],'svg')

close all
