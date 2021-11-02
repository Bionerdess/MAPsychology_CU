%% identify duration of excitation and inhibition periods

clear all

% data Directory

dataDir = 'C:\Users\LabAdmin\Dropbox\ANALYSIS\Anna-Lena\NAccVTA_April2021\'; % combined Data location
saveDir = 'C:\Users\LabAdmin\Dropbox\ANALYSIS\Anna-Lena\NAccVTA_April2021\images\';

statsCS = load([dataDir,'ExcInhPeriods_AB_CS.mat']);
statsCS = statsCS.CSstat;
statsUS = load([dataDir,'ExcInhPeriods_AB_US.mat']);
statsUS = statsUS.USstat;

ExcCS = nan(length(statsCS),2);
InhCS = nan(length(statsCS),2);
ExcUS = nan(length(statsUS),2);
InhUS = nan(length(statsUS),2);

for d = 1:length(statsCS)
    ExcCS(d,1) = sum([statsCS{1,d}]==1); % number of excited periods to A
    ExcCS(d,2) = sum([statsCS{2,d}]==1); % number of excited periods to B
    InhCS(d,1) = sum([statsCS{1,d}]==-1); % number of inhibited periods to A
    InhCS(d,2) = sum([statsCS{2,d}]==-1); % number of inhibited periods to B
    
    ExcUS(d,1) = sum([statsUS{1,d}]==1); % number of excited periods to A
    ExcUS(d,2) = sum([statsUS{2,d}]==1); % number of excited periods to B
    InhUS(d,1) = sum([statsUS{1,d}]==-1); % number of inhibited periods to A
    InhUS(d,2) = sum([statsUS{2,d}]==-1); % number of inhibited periods to B
end

ExcCS = ExcCS*100/1000; % transform number of 100ms bins into duration in seconds
InhCS = InhCS*100/1000; % transform number of 100ms bins into duration in seconds
ExcUS = ExcUS*100/1000; % transform number of 100ms bins into duration in seconds
InhUS = InhUS*100/1000; % transform number of 100ms bins into duration in seconds

fh1 = figure;
figure(fh1)
ecdf(ExcCS(:,1),'function','survivor')
hold on
ecdf(ExcCS(:,2),'function','survivor')
xlim([0 10])
ylabel('Proportion of neurons')
xlabel('Total duration of excitation [sec]')
legend('A','B','location','southeast')

fh2 = figure;
figure(fh2)
ecdf(InhCS(:,1),'function','survivor')
hold on
ecdf(InhCS(:,2),'function','survivor')
xlim([0 10])
ylabel('Proportion of neurons')
xlabel('Total duration of inhibition [sec]')
legend('A','B','location','southeast')

[hE,pE,statE] = kstest2(ExcCS(:,1),ExcCS(:,2));
[hI,pI,statI] = kstest2(InhCS(:,1),InhCS(:,2));

saveas(fh1,[saveDir,'ExcInh_DurationExc_survival.jpeg'],'jpeg')
saveas(fh1,[saveDir,'ExcInh_DurationExc_survival.svg'],'svg')
saveas(fh2,[saveDir,'ExcInh_DurationInh_survival.jpeg'],'jpeg')
saveas(fh2,[saveDir,'ExcInh_DurationInh_survival.svg'],'svg')

close all