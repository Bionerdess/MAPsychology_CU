%% Firing rate analysis - FR correlations CS & US for 4 500ms epochs

clear all

% get data dir
dataDir = 'C:\Users\LabAdmin\Dropbox\ANALYSIS\Anna-Lena\NAccVTA_April2021\'; % combined Data location
saveDir = 'C:\Users\LabAdmin\Dropbox\ANALYSIS\Anna-Lena\NAccVTA_April2021\images\';
groups = {'stage3'};

type = {'A_nr','B_nr'}; %

if smoothed ==1
    FRdata = load([dataDir,'FR_smoothed_stage3.mat']);
    CSresponse = load([dataDir,'CSresponse_smoothedDiffScore_2SD_stage3_500ms']);
else
    FRdata = load([dataDir,'FR_stage3.mat']);
end
NClass = load([dataDir,'VTA_values_2SD_',groups{1},'.mat']);
NClass = NClass.VTA_values;
Dind = find(cellfun(@(NClass) strcmp(NClass,'d'),NClass(:,5))==1);

uniqueID = cell2mat(FRdata.A.VTA.uniqueID);
count=1;

epoch=20;

Rvalues = nan(epoch,length(type));
Pvalues = nan(epoch,length(type));
for e=1:epoch
    avgCS1=nan(length(Dind),length(type));
    avgUS1 = nan(length(Dind),length(type));
    for f=1:length(type)
        if strcmp([type{f}],'A_nr')
            USoff=5;
        elseif strcmp([type{f}],'B_nr')
            USoff=0;
        end
        fh1 = figure;
        sig1=nan(length(Dind),4);
        for d=1:length(Dind)
            
            CS_z = FRdata.(type{f}).VTA.CS_z{Dind(d),1};
            US_z = FRdata.(type{f}).VTA.US_z{Dind(d),1};
            FR_CS1 = mean(CS_z(:,e*5-(5-1):e*5),1);
            FR_US1 = mean(US_z(:,USoff+1:USoff+5),1);
            
            avgCS1(d,f) = mean(FR_CS1);
            avgUS1(d,f) = mean(FR_US1);
            sig1(d,:)=0;
            [h1,p1,ci1,stats1] = ttest(FR_CS1,0);
            if p1<0.05 && stats1.tstat<0
                sig1(d,1)=1;
            end
            [h3,p3,ci3,stats3] = ttest(FR_US1,0);
            if p3<0.05
                sig1(d,3)=1;
            end
            
        end
        CS = rmmissing(avgCS1(:,f));
        U1 = rmmissing(avgUS1(:,f));
        sig1 = rmmissing(sig1);
        
        SD=1;
        IQR=0;
        if count==1
            if SD==1
                [~,zscore,idx1] = findOutliersSD(CS,1);
                [~,~,idx2] = findOutliersSD(U1,1);
                idx=[idx1;idx2];
            elseif IQR==1
                [~,idx1] = findOutliersIQR(CS,1);
                [~,idx2] = findOutliersIQR(U1,1);
                idx=[idx1;idx2];
            else
                idx=[];
            end
        end
        CS(idx,:)=[];
        U1(idx,:)=[];
        sig1(idx,:)=[];
        
        CS1(:,f) = CS;
        US1(:,f) = U1;
        
        [nCells,~]=size(CS1);
        
        for c=1:nCells
            figure(fh1)
            plot(CS1(c,f),US1(c,f),'.k','Markersize',25)
            hold on
            
        end
        count=count+1;
        
        sigCS1US1=find(sig1(:,1)==1);
        [rCS1US1,pCS1US1] = corr(CS1(:,f),US1(:,f),'Type','Spearman'); % CS1 vs US1 A trials - sigCS1US1
        Rvalues(e,f) = rCS1US1;
        Pvalues(e,f) = pCS1US1;
        
        y = [0 0 0];
        x = [-1 0 1];
        
        mdl1 = fitlm(CS1(:,f),US1(:,f));
        
        figure(fh1)
        legend(['r: ' num2str(rCS1US1)],['p: ' num2str(pCS1US1)])
        xlabel('CS')
        ylabel('US omit')
        title(['CS vs US (0-500ms) response ',type{f},' (',num2str((5*e-5)*100),'-',num2str((5*e)*100),' ms)'])
        
        saveas(fh1,[saveDir,'CorrelationSpearmanCS-US500_2SD_',type{f},'_',num2str((5*e-5)*100),'-',num2str((5*e)*100),'ms_SD_ALLneurons.svg'],'svg')
        saveas(fh1,[saveDir,'CorrelationSpearmanCS-US500_2SD_',type{f},'_',num2str((5*e-5)*100),'-',num2str((5*e)*100),'ms_SD_ALLneurons.jpeg'],'jpeg')
        
        
        
    end

    clear CS1
    
end



close all
