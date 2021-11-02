%% plot CDF of avg FR A&B 0-500ms

clear all

% get dataDir
dataDir = 'C:\Users\LabAdmin\Dropbox\ANALYSIS\Anna-Lena\NAccVTA_April2021\'; % combined Data location
saveDir = 'C:\Users\LabAdmin\Dropbox\ANALYSIS\Anna-Lena\NAccVTA_April2021\images\';

FR = load([dataDir,'FR_stage3']);
CellInfo = load([dataDir,'CellInfo_stage3']); 
NClass = load([dataDir,'VTA_values_2SD_stage3.mat']);
NClass = NClass.VTA_values;
Dind = find(cellfun(@(NClass) strcmp(NClass,'d'),NClass(:,5))==1);

type = {'A_nr','B_nr'};

count=1;

for f=1:length(type)
    if strcmp([type{f}],'A_nr')
        USoff=5;
    elseif strcmp([type{f}],'B_nr')
        USoff=0;
    end
    
    avgCS = nan(length(Dind),1);
    avgUS = nan(length(Dind),1);
    
    for c = 1:length(Dind)
        CS = FR.(type{f}).VTA.CS_z{Dind(c),1};
        US = FR.(type{f}).VTA.US_z{Dind(c),1};
        
        CS = mean(CS(:,1:10),2);
        avgCS(c) = mean(CS);
        US = mean(US(:,USoff+1:USoff+5),2);
        avgUS(c) = mean(US);
    end
    SD=1;
    IQR=0;
    if count==1
        if SD==1
            [~,zscore,idx1] = findOutliersSD(avgCS,1);
            [~,~,idx2] = findOutliersSD(avgUS,1);
            idx=[idx1;idx2];
        elseif IQR==1
            [~,idx1] = findOutliersIQR(avgCS,1);
            [~,idx2] = findOutliersIQR(avgUS,1);
            idx=[idx1;idx2];
        else
            idx=[];
        end
    end
    avgCS(idx,:)=[];
    avgUS(idx,:)=[];
    
    C(:,f) = avgCS;
    U(:,f) = avgUS;
    
    count=count+1;
end

fh1 = figure;
figure(fh1)
ecdf(C(:,1))
hold on
ecdf(C(:,2))
legend('A','B')
title('CDF - avgFR CS (0-500ms)')
ylabel('Proportion of Neurons')
xlabel('avgFR [Hz]')

saveas(fh1,[saveDir,'CDF_avgFRzscore_CS_0-1sec.svg'],'svg')
[h,p,stats] = kstest2(C(:,1),C(:,2));
disp('A vs. B')
p
stats

fh2 = figure;
figure(fh2)
ecdf(U(:,1))
hold on
ecdf(U(:,2))
legend('A-US','B-US')
title('CDF - avgFR US omission (0-500ms)')
ylabel('Proportion of Neurons')
xlabel('avgFR [Hz]')

saveas(fh2,[saveDir,'CDF_avgFR_US_0-500ms.svg'],'svg')
[h,p,stats] = kstest2(U(:,1),U(:,2));
disp('A US vs. B US')
p
stats