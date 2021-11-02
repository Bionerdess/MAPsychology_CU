%% Generate Exc/Inh stats CS & US

clear all

% data Directory

dataDir = 'C:\Users\LabAdmin\Dropbox\ANALYSIS\Anna-Lena\NAccVTA_April2021\'; % combined Data location
saveDir = 'C:\Users\LabAdmin\Dropbox\ANALYSIS\Anna-Lena\NAccVTA_April2021\images\';

FR = load([dataDir,'FR_stage3']);
CellInfo = load([dataDir,'CellInfo_stage3']); % load CellInfo data (animal,CellID,path,peakCh,session,stage)
NClass = load([dataDir,'VTA_values_2SD_stage3.mat']);
NClass = NClass.VTA_values;
Dind = find(cellfun(@(NClass) strcmp(NClass,'d'),NClass(:,5))==1);

type = {'A_nr','B_nr'};

CSstat = cell(length(type),length(Dind));
USstat = cell(length(type),length(Dind));
for f = 1:length(type)
    if strcmp([type{f}],'A_nr')
        USoff=5;
    elseif strcmp([type{f}],'B_nr')
        USoff=0;
    end
    for d = 1:length(Dind)
        BL = FR.(type{f}).VTA.baseline{Dind(d),1};
        CS = FR.(type{f}).VTA.CS{Dind(d),1};
        US = FR.(type{f}).VTA.US{Dind(d),1};
        
        [statsCS,cumSumCS,avgCumBL] = calcCumSum(CS,BL,'alpha',0.001,'tail',2);
        CSstat{f,d} = statsCS;
        
        [statsUS,cumSumUS,avgCumBL2] = calcCumSum(US,BL(:,1:size(US,2)),'alpha',0.001,'tail',2);
        USstat{f,d} = statsUS;
    end
end
save([dataDir,'ExcInhPeriods_AB_CS'],'CSstat');
save([dataDir,'ExcInhPeriods_AB_US'],'USstat');

%% Identify periods of excitation and inhibition during 0-3sec of CS

clear all

% data Directory

dataDir = 'C:\Users\LabAdmin\Dropbox\ANALYSIS\Anna-Lena\NAccVTA_April2021\'; % combined Data location
saveDir = 'C:\Users\LabAdmin\Dropbox\ANALYSIS\Anna-Lena\NAccVTA_April2021\images\';

FR = load([dataDir,'FR_stage3']);
CellInfo = load([dataDir,'CellInfo_stage3']); % load CellInfo data (animal,CellID,path,peakCh,session,stage)
NClass = load([dataDir,'VTA_values_2SD_stage3.mat']);
NClass = NClass.VTA_values;
Dind = find(cellfun(@(NClass) strcmp(NClass,'d'),NClass(:,5))==1);

statsCS = load([dataDir,'ExcInhPeriods_AB_CS.mat']);
statsCS = statsCS.CSstat;
statsUS = load([dataDir,'ExcInhPeriods_AB_US.mat']);
statsUS = statsUS.USstat;

stim = {'CS','US'};

window = 20; % set number of 100ms bins considered
nSig = 3; % set number of bins that have to be significantly inhibited (-1) or excited (+1)

statsAll = {statsCS}; % ,statsUS

for s = 1:length(statsAll)
    if s==2 % US
        start = 5; % US off for A trials is at 10.5 sec -> 5 bins into US
        window = 43;
    else
        start=1;
    end
    stats = statsAll{s};
    res = cell(2,size(stats,2));
    Binh = nan(length(stats),2);
    Aexc = nan(length(stats),2);
    Bexc = nan(length(stats),2);
    Ainh = nan(length(stats),2);
    
    for c = 1:length(stats)
        uniqueID = NClass{Dind(c),3};
        Binh(c,2) = uniqueID;
        Aexc(c,2) = uniqueID;
        Bexc(c,2) = uniqueID;
        Ainh(c,2) = uniqueID;
        Asum = nan(1,window);
        Bsum = nan(1,window);
        
        for ia = start:window+start
            Asum(1,ia-start+1) = sum(stats{1,c}(ia:ia+nSig-1));
        end
        for i = 1:window
            Bsum(1,i) = sum(stats{2,c}(i:i+nSig-1));
        end
        res{1,c}{1} = find(Asum ==nSig);
        res{1,c}{2} = find(Asum ==-nSig);
        res{2,c}{1} = find(Bsum ==nSig);
        res{2,c}{2} = find(Bsum ==-nSig);
        
        if ~isempty(res{2,c}{2})
            Binh(c,1) = res{2,c}{2}(1);
            % uncomment for "out of all neurons"
%         else
%             Binh(c,1) = 100;
        end
        if ~isempty(res{1,c}{1})
            Aexc(c,1) = res{1,c}{1}(1);
%         else
%             Aexc(c,1) = 100;
        end
        if ~isempty(res{2,c}{1})
            Bexc(c,1) = res{2,c}{1}(1);
%         else
%             Bexc(c,1) = 100;
        end
        if ~isempty(res{1,c}{2})
            Ainh(c,1) = res{1,c}{2}(1);
%         else
%             Ainh(c,1) = 100;
        end
    end
    
    
    indexcAexcB = find(~isnan(Aexc(:,1))==1 & ~isnan(Bexc(:,1)) ==1 & ~isnan(Ainh(:,1))==0 & ~isnan(Binh(:,1))==0);
    indexcAinhB = find(~isnan(Aexc(:,1))==1 & ~isnan(Bexc(:,1)) ==0 & ~isnan(Ainh(:,1))==0 & ~isnan(Binh(:,1))==1);
    indinhAexcB = find(~isnan(Aexc(:,1))==0 & ~isnan(Bexc(:,1)) ==1 & ~isnan(Ainh(:,1))==1 & ~isnan(Binh(:,1))==0);
    indinhAinhB = find(~isnan(Aexc(:,1))==0 & ~isnan(Bexc(:,1)) ==0 & ~isnan(Ainh(:,1))==1 & ~isnan(Binh(:,1))==1);
    
    cMat(1,1) = length(indexcAexcB)/length(stats);
    cMat(1,2) = length(indexcAinhB)/length(stats);
    cMat(2,1) = length(indinhAexcB)/length(stats);
    cMat(2,2) = length(indinhAinhB)/length(stats);
    
    fh1 = figure;
    figure(fh1)
    imagesc(cMat)
    colormap(jet)
    colorbar
    
    saveas(fh1,[saveDir,'CDF_ExcA-InhB_startBins1sec.svg'],'svg')

    
    fh3 = figure;
    figure(fh3)
    cdfplot(Aexc(:,1))
    hold on
    cdfplot(Ainh(:,1))
    legend('A exc','A inh','location','southeast')
    title(['CDF ',stim{s},' - Start bin of inhibition (A) or excitation (A) - 2 bins'])
    ylabel('Proportion of Neurons')
    xlabel('Start bin [sec]')
%     xlim([0 window])
    ylim([0 1])
    
    saveas(fh3,[saveDir,'CDF_ExcA-InhA_startBins',num2str(window),'bins2_',stim{s},'.svg'],'svg')
    [h,p,stats] = kstest2(Aexc(:,1),Ainh(:,1));
    disp('A exc vs. A inh')
    p
    stats
    
    fh4 = figure;
    figure(fh4)
    cdfplot(Bexc(:,1))
    hold on
    cdfplot(Binh(:,1))
    legend('B exc','B inh','location','southeast')
    title(['CDF ',stim{s},' - Start bin of inhibition (B) or excitation (B) - 2 bins'])
    ylabel('Proportion of Neurons')
    xlabel('Start bin [sec]')
%     xlim([0 window])
    ylim([0 1])
    
    saveas(fh4,[saveDir,'CDF_ExcB-InhB_startBins',num2str(window),'bins2_',stim{s},'.svg'],'svg')
    [h,p,stats] = kstest2(Bexc(:,1),Binh(:,1));
    disp('B exc vs. B inh')
    p
    stats
    
    fh5 = figure;
    figure(fh5)
    ecdf(Aexc(:,1),'Function','survivor')
    hold on
    ecdf(Bexc(:,1),'Function','survivor')
    legend('A exc','B exc','location','northeast')
    title(['CDF ',stim{s},' - Start bin of excitation (A) or excitation (B) - ',num2str(nSig),' bins'])
    ylabel('Proportion of Neurons')
    xlabel('Start bin [sec]')
%     xlim([0 window])
    ylim([0 1])
    
    saveas(fh5,[saveDir,'CDF_ExcA-ExcB_startBins_survival',num2str(window+nSig),'bins',num2str(nSig),'_',stim{s},'.svg'],'svg')
    [h,p,stats] = kstest2(Aexc(:,1),Bexc(:,1));
    disp('A exc vs. B exc')
    p
    stats
    
    fh6 = figure;
    figure(fh6)
    ecdf(Ainh(:,1),'Function','survivor')
    hold on
    ecdf(Binh(:,1),'Function','survivor')
    legend('A inh','B inh','location','northeast')
    title(['CDF ',stim{s},' - Start bin of inhibition (A) or inhibition (B) - ',num2str(nSig),' bins'])
    ylabel('Proportion of Neurons')
    xlabel('Start bin [sec]')
%     xlim([0 window])
    ylim([0 1])
    
    saveas(fh6,[saveDir,'CDF_InhA-InhB_startBins_survival',num2str(window+nSig),'bins',num2str(nSig),'_',stim{s},'.svg'],'svg')
    [h,p,stats] = kstest2(Ainh(:,1),Binh(:,1));
    disp('A inh vs. B inh')
    p
    stats
    
    close all
end