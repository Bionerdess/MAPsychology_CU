%% Plot PSTH A & B mean +/- SEM

clear all

% get data dir
dataDir = 'C:\Users\LabAdmin\Dropbox\ANALYSIS\Anna-Lena\NAccVTA_April2021\'; % combined Data location
saveDir = 'C:\Users\LabAdmin\Dropbox\ANALYSIS\Anna-Lena\NAccVTA_April2021\images\';

FR = load([dataDir,'FR_stage3']);
NClass = load([dataDir,'VTA_values_2SD_stage3']);
NClass = NClass.VTA_values;
Dind = find(cellfun(@(NClass) strcmp(NClass,'d'),NClass(:,5))==1);

type = {'A_nr','B_nr'};

acc = load([dataDir,'Accuracy_byNeuron_USdecoder']);
accuracy = acc.accuracy;
accuracyBL = acc.accuracyBL;

data = cell(1,length(type));
dataUS = cell(1,length(type));
for f=1:length(type)
    if strcmp(type{f},'A_nr')
        USoff = 6;
    elseif strcmp(type{f},'B_nr')
        USoff = 1;
    end
    for d = 1:length(Dind)
        baseline = FR.(type{f}).VTA.baseline_z{Dind(d),1};
        CS_z = FR.(type{f}).VTA.CS_z{Dind(d),1};
        US_z = FR.(type{f}).VTA.US_z{Dind(d),1};
        
        temp = [baseline(:,size(baseline,2)/2+1:size(baseline,2)),CS_z];
        data{f} = cat(1,data{f},mean(temp,1));
        dataUS{f} = cat(1,dataUS{f},mean(US_z(:,USoff:USoff+44),1));
    end
end

x1 = size(baseline,2)/2;
c = [0 128 102; 128 0 51]/255;
AB = cat(1,data{1},data{2});
ABUS = cat(1,dataUS{1},dataUS{2});
avgAB_CS = (data{1}+data{2})/2;
diffA_CS = data{1}-mean(AB,1);
diffB_CS = data{2}-mean(AB,1);
avgAB_US = (dataUS{1}+dataUS{2})/2;
diffA_US = dataUS{1}-mean(ABUS,1);
diffB_US = dataUS{2}-mean(ABUS,1);
diffCS{1} = diffA_CS;
diffCS{2} = diffB_CS;
diffUS{1} = diffA_US;
diffUS{2} = diffB_US;

[pop1,~] = find(sum(accuracy(:,3:8)<accuracyBL(:,3:8),2)>3); % population that decodes worse than random
[pop2,~] = find(sum(accuracy(:,3:8)>accuracyBL(:,3:8),2)>3); % population that decodes better than random
[pop3,~] = find(mean(accuracy(:,3:7),2)<mean(accuracyBL(:,3:7),2)& mean(accuracy(:,10:20),2)>mean(accuracyBL(:,10:20),2)); % population that switches decoding than random
[pop4,~] = find((mean(accuracy(:,3:7),2)>mean(accuracyBL(:,3:7),2)& mean(accuracy(:,10:20),2)>mean(accuracyBL(:,10:20),2)) | (mean(accuracy(:,3:7),2)<mean(accuracyBL(:,3:7),2)& mean(accuracy(:,10:20),2)<mean(accuracyBL(:,10:20),2))); % population that does not switch
pop1 = unique(pop1);
pop2 = unique(pop2);

datPop1{1} = data{1}(pop1,:);
datPop1{2} = data{2}(pop1,:);
datPop2{1} = data{1}(pop2,:);
datPop2{2} = data{2}(pop2,:);
datPop3{1} = data{1}(pop3,:);
datPop3{2} = data{2}(pop3,:);
datPop4{1} = data{1}(pop4,:);
datPop4{2} = data{2}(pop4,:);

datUS1{1} = dataUS{1}(pop1,:);
datUS1{2} = dataUS{2}(pop1,:);
datUS2{1} = dataUS{1}(pop2,:);
datUS2{2} = dataUS{2}(pop2,:);
datUS3{1} = dataUS{1}(pop3,:);
datUS3{2} = dataUS{2}(pop3,:);
datUS4{1} = dataUS{1}(pop4,:);
datUS4{2} = dataUS{2}(pop4,:);

mPop1(:,1) = mean(datPop1{1}(:,53:57),2);
mPop1(:,2) = mean(datPop1{2}(:,53:57),2);
mPop2(:,1) = mean(datPop2{1}(:,53:57),2);
mPop2(:,2) = mean(datPop2{2}(:,53:57),2);
mPop1US(:,1) = mean(datUS1{1}(:,1:5),2);
mPop1US(:,2) = mean(datUS1{2}(:,1:5),2);
mPop2US(:,1) = mean(datUS2{1}(:,1:5),2);
mPop2US(:,2) = mean(datUS2{2}(:,1:5),2);
mBL1(:,1) = mean(datPop1{1}(:,1:50),2);
mBL1(:,2) = mean(datPop1{2}(:,1:50),2);
mBL1 = mean(mBL1,2);
mBL2(:,1) = mean(datPop2{1}(:,1:50),2);
mBL2(:,2) = mean(datPop2{2}(:,1:50),2);
mBL2 = mean(mBL2,2);


mPop1 = (mPop1(:,1)-mPop1(:,2))/2;
mPop1US = (mPop1US(:,1)-mPop1US(:,2))/2;




[ra1,pa1] = corr(mPop1(:,1),mPop1US(:,1),'Type','Spearman');
[ra2,pa2] = corr(mPop2(:,1),mPop2US(:,1),'Type','Spearman');

mdl1 = fitlm(mPop1(:,1),mPop1US(:,1));
mdl3 = fitlm(mPop2(:,1),mPop2US(:,1));

fh1 = figure;
figure(fh1)
plot(mPop1(:,1),mPop1US(:,1),'k.','MarkerSize',30)
hold on
% plot(mdl1)
legend(['p: ',num2str(pa1),'; r: ',num2str(ra1)])
xlabel('CS')
ylabel('US omission')
title('Pop Worse - AB (inverse avg)')

saveas(fh1,[saveDir,'FINAL_Correlation_popWorse4_AB3-8.svg'],'svg')


fh3 = figure;
figure(fh3)
plot(mPop2(:,1),mPop2US(:,1),'k.','MarkerSize',30)
hold on
% plot(mdl3)
legend(['p: ',num2str(pa2),'; r: ',num2str(ra2)])
xlabel('CS')
ylabel('US omission')
title('Pop Better - AB (inverse avg)')

saveas(fh3,[saveDir,'FINAL_Correlation_popBetter4_AB3-8.svg'],'svg')

mA1 = mean(datPop1{1}(:,53:58),2);
mB1 = mean(datPop1{2}(:,53:58),2);
mA2 = mean(datPop2{1}(:,53:58),2);
mB2 = mean(datPop2{2}(:,53:58),2);
mUSA1 = mean(datUS1{1}(:,1:5),2);
mUSB1 = mean(datUS1{2}(:,1:5),2);
mUSA2 = mean(datUS2{1}(:,1:5),2);
mUSB2 = mean(datUS2{2}(:,1:5),2);

[pABw,~,statABw] = ranksum(mA1,mB1); % non-parametric Mann-Whitney-U test
[~,pABw2,ciABw2,statABw2] = ttest2(mA1,mB1);
[pAw,~,statAw] = ranksum(mA1,mBL1);
[~,pAw2,ciAw2,statAw2] = ttest2(mA1,mBL1);
[pBw,~,statBw] = ranksum(mB1,mBL1);
[~,pBw2,ciBw2,statBw2] = ttest2(mB1,mBL1);
median(mA1)
median(mB1)
median(mBL1)

[pABb,~,statABb] = ranksum(mA2,mB2); % non-parametric Mann-Whitney-U test
[~,pABb2,ciABb2,statABb2] = ttest2(mA2,mB2);
[pAb,~,statAb] = ranksum(mA2,mBL2);
[~,pAb2,ciAb2,statAb2] = ttest2(mA2,mBL2);
[pBb,~,statBb] = ranksum(mB2,mBL2);
[~,pBb2,ciBb2,statBb2] = ttest2(mB2,mBL2);
median(mA2)
median(mB2)
median(mBL2)

[pABUw,~,statABUw] = ranksum(mUSA1,mUSB1); % non-parametric Mann-Whitney-U test
[~,pABUw2,ciABUw2,statABUw2] = ttest2(mUSA1,mUSB1);
[pUAw,~,statUAw] = ranksum(mUSA1,mBL1);
[~,pAUw2,ciAUw2,statAUw2] = ttest2(mUSA1,mBL1);
[pBUw,~,statBUw] = ranksum(mUSB1,mBL1);
[~,pBUw2,ciBUw2,statBUw2] = ttest2(mUSB1,mBL1);
median(mUSA1)
median(mUSB1)

[pABUb,~,statABUb] = ranksum(mUSA2,mUSB2); % non-parametric Mann-Whitney-U test
[~,pABUb2,ciABUb2,statABUb2] = ttest2(mUSA2,mUSB2);
[pAUb,~,statAUb] = ranksum(mUSA2,mBL2);
[~,pAUb2,ciAUb2,statAUb2] = ttest2(mUSA2,mBL2);
[pBUb,~,statBUb] = ranksum(mUSB2,mBL2);
[~,pBUb2,ciBUb2,statBUb2] = ttest2(mUSB2,mBL2);
median(mUSA2)
median(mUSB2)


[X,fig1] = PSTH(datPop1,'vertline',[x1 x1+100],'color',c);
figure(fig1)
ylim([-2 12])
ylabel('FR [z-score]')
xlabel('Time [ms]')
title('Worse than random model')
fig1.Position = [50 50 800 400];
saveas(fig1,[saveDir,'FINAL_PSTH_WorseGLM3-8_2SD_AB_BLCS34.svg'],'svg')

[X,fig2] = PSTH(datPop2,'vertline',[x1 x1+100],'color',c);
figure(fig2)
ylim([-2 12])
ylabel('FR [z-score]')
xlabel('Time [ms]')
title('Better than random model')
fig2.Position = [50 50 800 400];
saveas(fig2,[saveDir,'FINAL_PSTH_BetterGLM3-8_2SD_AB_BLCS4.svg'],'svg')

[XUS,figUS1] = PSTH(datUS1,'vertline',[x1 x1+100],'color',c);
figure(figUS1)
ylabel('FR [z-score]')
xlabel('Time [ms]')
ylim([-2 12])
title('Worse than random model')
figUS1.Position = [50 50 800 400];
saveas(figUS1,[saveDir,'FINAL_PSTH_WorseGLM_2SD_AB_US4.svg'],'svg')

[XUS,figUS2] = PSTH(datUS2,'vertline',[x1 x1+100],'color',c);
figure(figUS2)
ylabel('FR [z-score]')
xlabel('Time [ms]')
ylim([-2 12])
title('Better than random model')
figUS2.Position = [50 50 800 400];
saveas(figUS2,[saveDir,'FINAL_PSTH_BetterGLM_2SD_AB_US4.svg'],'svg')

close all