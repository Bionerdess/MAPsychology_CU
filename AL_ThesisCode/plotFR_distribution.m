%% Distribution of FR 0-500ms CS  & US (from US off)

clear all

% get data dir
dataDir = 'C:\Users\LabAdmin\Dropbox\ANALYSIS\Anna-Lena\NAccVTA_April2021\'; % combined Data location
saveDir = 'C:\Users\LabAdmin\Dropbox\ANALYSIS\Anna-Lena\NAccVTA_April2021\';

FRdata = load([dataDir,'FR_','stage3','.mat']);

NClass = load([dataDir,'VTA_values_2SD_','stage3','.mat']);
NClass = NClass.VTA_values;
Dind = find(cellfun(@(NClass) strcmp(NClass,'d'),NClass(:,5))==1);
type = {'A_nr','B_nr'}; 

epoch=4;
count=1;
p_CS = nan(epoch,length(type));
t_CS = nan(epoch,length(type));
p_US = nan(epoch,length(type));
t_US = nan(epoch,length(type));
CS = cell(length(type),1);
U1 = cell(length(type),1);
for e=1:epoch
    dfR = nan(length(Dind),length(type)*2);
    tR = nan(length(Dind),length(type)*2);
    pR = nan(length(Dind),length(type)*2);
    
    avgFRCS1 = nan(length(Dind),length(type)*2);
    avgFRUS1 = nan(length(Dind),length(type)*2);
    avgFRCS2 = nan(length(Dind),length(type)*2);
    avgFRUS2 = nan(length(Dind),length(type)*2);
    for f = 1:length(type)
        if strcmp([type{f}],'A_nr')
            USoff=5;
        elseif strcmp([type{f}],'B_nr')
            USoff=0;
        end
        for d=1:length(Dind)
                CS_z = FRdata.(type{f}).VTA.CS_z{Dind(d),1};
                US_z = FRdata.(type{f}).VTA.US_z{Dind(d),1};
                FRCS1 = mean(CS_z(:,e*5-(5-1):e*5),1); % 5 bins=500ms
                FRUS1 = mean(US_z(:,USoff+1:USoff+5),1);
                avgFRCS1(d,f*2-1) = mean(FRCS1);
                avgFRUS1(d,f*2-1) = mean(FRUS1);
                
                [hc1,pc1,cic1,statsc1] = ttest(FRCS1,0);
                if pc1<0.05 && statsc1.tstat <0
                    avgFRCS1(d,f*2)=1;
                elseif pc1<0.05 && statsc1.tstat >0
                    avgFRCS1(d,f*2)=2;
                else
                    avgFRCS1(d,f*2)=0;
                end
                [hu1,pu1,ciu1,statsu1] = ttest(FRUS1,0);
                if pu1<0.05 && statsu1.tstat <0
                    avgFRUS1(d,f*2)=1;
                elseif pu1<0.05 && statsu1.tstat >0
                    avgFRUS1(d,f*2)=2;
                else
                    avgFRUS1(d,f*2)=0;
                end
                
                pR(d,f*2-1) = pc1;
                tR(d,f*2-1) = statsc1.tstat;
                dfR(d,f*2-1) = statsc1.df;
                pR(d,f*2) = pu1;
                tR(d,f*2) = statsu1.tstat;
                dfR(d,f*2) = statsu1.df;
        end
        
        CS{f} = rmmissing(avgFRCS1(:,f*2-1:f*2));
        U1{f} = rmmissing(avgFRUS1(:,f*2-1:f*2));
        
       SD=1;
        IQR=0;
        if count==1
            if SD==1
                [~,zscore,idx1] = findOutliersSD(CS{f},1);
                [~,~,idx2] = findOutliersSD(U1{f},1);
                idx=[idx1;idx2];
            elseif IQR==1
                [~,idx1] = findOutliersIQR(CS{f},1);
                [~,idx2] = findOutliersIQR(U1{f},1);
                idx=[idx1;idx2];
            else
                idx=[];
            end
        end
        CS{f}(idx,:)=[];
        U1{f}(idx,:)=[];
        
        CS1(:,1:2) = CS{f};
        US1(:,1:2) = U1{f};
        excCS1{:,f} = find(CS1(:,2)==2);
        inhCS1{:,f} = find(CS1(:,2)==1);
        nsCS1{:,f} = find(CS1(:,2)==0);
        excUS1{:,f} = find(US1(:,2)==2);
        inhUS1{:,f} = find(US1(:,2)==1);
        nsUS1{:,f} = find(US1(:,2)==0);
        
        [hCS,pCS,ciCS,statsCS] = ttest(CS1(:,1),0);
        [hUS,pUS,ciUS,statsUS] = ttest(US1(:,1),0);
        
        p_CS(e,f) = pCS;
        t_CS(e,f) = statsCS.tstat;
        p_US(e,f) = pUS;
        t_US(e,f) = statsUS.tstat;
        
        fh1 = figure;
        fh2=figure;
        
        figure(fh1)
        histogram(CS1(:,1),'BinWidth',range(CS1(:,1))/15,'FaceColor','k')
        hold on
        plot([0 0],[0 15],'k-')
        xlabel('FR [zScore]')
        ylabel('# DA neurons')
        title(['Population avgFR - ',type{f},' CS (',num2str((5*e-5)*100),'-',num2str((5*e)*100),' ms)'])
        
        figure(fh2)
        histogram(US1(:,1),'BinWidth',range(US1(:,1))/15,'FaceColor','k')
        hold on
        plot([0 0],[0 15],'k-')
        xlabel('FR [zScore]')
        ylabel('# DA neurons')
        title(['Population avgFR - ',type{f},' US (',num2str((5*e-5)*100),'-',num2str((5*e)*100),' ms)'])
        
        saveas(fh1,[saveDir,'images\DistributionPopulation_avgFR_ALLneurons_SD_',type{f},'_CS',num2str((5*e-5)*100),'-',num2str((5*e)*100),'ms.svg'],'svg')
        saveas(fh2,[saveDir,'images\DistributionPopulation_avgFR__ALLneurons_SD_',type{f},'_US',num2str((5*e-5)*100),'-',num2str((5*e)*100),'ms.svg'],'svg')
        count=count+1;
    end

    close all
end