%% CS decoder based on top PCs, compared to LDA
% This protocol runs a PCA for each epoch of CS, then runs an LDA and GLM
% on the same epoch.

clear all

% get data dir
dataDir = 'C:\Users\LabAdmin\Dropbox\ANALYSIS\Anna-Lena\NAccVTA_April2021\'; % combined Data location
saveDir = 'C:\Users\LabAdmin\Dropbox\ANALYSIS\Anna-Lena\NAccVTA_April2021\images\';

FRdata = load([dataDir,'FR_stage3']);
NClass = load([dataDir,'VTA_values_2SD_stage3.mat']);
NClass = NClass.VTA_values;
Dind = find(cellfun(@(NClass) strcmp(NClass,'d'),NClass(:,5))==1);

type = {'A_nr','B_nr'};

select = 40;
reps = 1000; 
smpl = randi(length(Dind),reps,select);

epoch = 100;
AUCg = nan(reps,epoch);
AUCl = nan(reps,epoch);
AUCgShuff = nan(reps,epoch);
AUClShuff = nan(reps,epoch);

GLMacc = nan(reps,epoch,2);
LDAacc = nan(reps,epoch,2);
GLMaccShuff = nan(reps,epoch,2);
LDAaccShuff = nan(reps,epoch,2);

accuracy = struct;
AUC = struct;
VarExp = nan(15,epoch,reps);
PCcoeff = cell(4,1);

tic
for comb = 1:reps
    Y = cat(1,zeros(8,1),ones(8,1)); % trial types (8xA: 0; 8xB: 1)
    s = randperm(16); % shuffle trials, keep the same for all epochs
    Y= Y(s,:); % shuffle Y the same as trial FR data

    c = cvpartition(Y,'HoldOut',0.5); % split train/test data the same for all epochs
    idxTrain = training(c,1);
    idxTest = ~idxTrain;
    
    for e=1:epoch
        neuronData = [];
        BLdata = [];
        for f = 1:length(type)
            temp = nan(8,select);
            tempBL = nan(8,select);
            for d = 1:select
                CS = FRdata.(type{f}).VTA.CS_z{Dind(smpl(comb,d)),1};
                BL = FRdata.(type{f}).VTA.baseline_z{Dind(smpl(comb,d)),1};
                tempBL(:,d) = BL(:,e);
                temp(:,d) = CS(:,e);
            end
            
            neuronData = cat(1,neuronData,temp);
            BLdata = cat(1,BLdata,tempBL);
        end
        % shuffle trials
        neuronData = neuronData(s,:);
        centeredData = neuronData-mean(neuronData,1);
        BLdata = BLdata(s,:);
        centeredBL = BLdata-mean(BLdata,1);
        
        % run PCA
        [coeff,score,~,~,varExp] = pca(neuronData); %
        VarExp(:,e,comb) = varExp;
        PCcoeff{1}(:,e,comb) = coeff(:,1);
        PCcoeff{2}(:,e,comb) = coeff(:,2);
        PCcoeff{3}(:,e,comb) = coeff(:,3);
        PCcoeff{4}(:,e,comb) = coeff(:,4);
        % dot product of top 8 PCs and mean centered data
        dotDat = centeredData*coeff(:,1:4);
        dotBL = centeredBL*coeff(:,1:4);
        
        % split into train and test data
        XTrain = dotDat(idxTrain,:);
        yTrain = Y(idxTrain);
        XTest = dotDat(idxTest,:);
        yTest = Y(idxTest);
        
        
        % train LDA 
        mdl = fitcdiscr(XTrain,yTrain);
        % test LDA on test trials
        pred = predict(mdl,XTest);
        cm = confusionmat(yTest,pred);
        [~,~,~,AUCl(comb,e)] = perfcurve(yTest,pred,1);
        
        LDAacc(comb,e,1) = cm(1,1)/(cm(1,1)+cm(1,2)); % percent accurate out of A trials
        LDAacc(comb,e,2) = cm(2,2)/(cm(2,1)+cm(2,2));
        
        % test validation LDA on BL trials
        predShuff = predict(mdl,dotBL);
        cs = confusionmat(Y,predShuff);
        [~,~,~,AUClShuff(comb,e)] = perfcurve(Y,predShuff,1);
        
        LDAaccShuff(comb,e,1) = cs(1,1)/(cs(1,1)+cs(1,2)); % percent accurate out of A trials
        LDAaccShuff(comb,e,2) = cs(2,2)/(cs(2,1)+cs(2,2));
        
        % train LASSO GLM on top 8 PCA scores
        [B,FitInfo] = lassoglm(XTrain,yTrain,'binomial','CV',8,'Alpha',.5); 
        idxLambdaMinDeviance = FitInfo.IndexMinDeviance;
        B0 = FitInfo.Intercept(idxLambdaMinDeviance);
        coef = [B0; B(:,idxLambdaMinDeviance)];
        
        % test GLM 
        yhat = glmval(coef,XTest,'logit');
        yHatBinom = double(yhat>=0.5);
        
        [X,Y1,T,AUCg(comb,e)] = perfcurve(yTest,yhat,1);
        
        mat = confusionmat(yTest,yHatBinom);
        
        GLMacc(comb,e,1) = mat(1,1)/(mat(1,1)+mat(1,2));
        GLMacc(comb,e,2) = mat(2,2)/(mat(2,1)+mat(2,2));
        
        % test validation GLM on shuffled test data 
        yhatShuff = glmval(coef,dotBL,'logit');
        yHatBinomShuff = double(yhatShuff>=0.5);
        
        [~,~,~,AUCgShuff(comb,e)] = perfcurve(Y,yhatShuff,1);
        
        matShuff = confusionmat(Y,yHatBinomShuff);
        
        GLMaccShuff(comb,e,1) = matShuff(1,1)/(matShuff(1,1)+matShuff(1,2));
        GLMaccShuff(comb,e,2) = matShuff(2,2)/(matShuff(2,1)+matShuff(2,2));
        
    end
    toc
end

mGLMacc = reshape(mean(GLMacc,1),size(GLMacc,2),size(GLMacc,3));
mLDAacc = reshape(mean(LDAacc,1),size(LDAacc,2),size(LDAacc,3));
mGLMaccShuff = reshape(mean(GLMaccShuff,1),size(GLMaccShuff,2),size(GLMaccShuff,3));
mLDAaccShuff = reshape(mean(LDAaccShuff,1),size(LDAaccShuff,2),size(LDAaccShuff,3));

fh1 = figure;
figure(fh1)
plot(mean(mGLMacc,2),'LineWidth',2)
hold on
plot(mean(mGLMaccShuff,2),'LineWidth',2)
hold on
plot([0 epoch],[.5 .5],'k--')
ylim([0 1])
legend('GLMacc','validation GLMacc')
xlabel('Time (100ms bins)')
ylabel('avg % accurate')

fh2 = figure;
figure(fh2)
plot(mean(mLDAacc,2),'LineWidth',2)
hold on
plot(mean(mLDAaccShuff,2),'LineWidth',2)
hold on
plot([0 epoch],[.5 .5],'k--')
ylim([0 1])
legend('LDAacc','validation LDAacc')
xlabel('Time (100ms bins)')
ylabel('avg % accurate')

fh3 = figure;
figure(fh3)
plot(mean(AUCg,1),'LineWidth',3)
hold on
plot(mean(AUCl,1),'Linewidth',3)
hold on
plot([0 epoch],[.5 .5],'k--')
ylim([0 1])
legend('GLM','LDA')
xlabel('Time (100ms bins)')
ylabel('AUC')

fh4 = figure;
figure(fh4)
plot(mean(mGLMacc,2),'LineWidth',2)
hold on
plot(mean(mLDAacc,2),'LineWidth',2)
hold on
plot([0 epoch],[.5 .5],'k--')
ylim([0 1])
legend('GLM','LDA')
xlabel('Time (100ms bins)')
ylabel('avg % accurate')

fh5 = figure;
figure(fh5)
plot(mean(AUCg,1),'LineWidth',3)
hold on
plot(mean(AUCgShuff,1),'Linewidth',3)
hold on
plot([0 epoch],[.5 .5],'k--')
ylim([0 1])
legend('GLM','validation GLM')
xlabel('Time (100ms bins)')
ylabel('AUC')

fh6 = figure;
figure(fh6)
plot(mean(AUCl,1),'LineWidth',3)
hold on
plot(mean(AUClShuff,1),'Linewidth',3)
hold on
plot([0 epoch],[.5 .5],'k--')
ylim([0 1])
legend('LDA','validation LDA')
xlabel('Time (100ms bins)')
ylabel('AUC')


saveas(fh1,[saveDir,'avgAccuracyGLM_CSdecoder_PCA_',num2str(epoch),'epochs_validatedZ4.jpeg'],'jpeg')
saveas(fh1,[saveDir,'avgAccuracyGLM_CSdecoder_PCA_',num2str(epoch),'epochs_validatedZ4.svg'],'svg')
saveas(fh2,[saveDir,'avgAccuracyLDA_CSdecoder_PCA_',num2str(epoch),'epochs_validatedZ4.jpeg'],'jpeg')
saveas(fh2,[saveDir,'avgAccuracyLDA_CSdecoder_PCA_',num2str(epoch),'epochs_validatedZ4.svg'],'svg')
saveas(fh3,[saveDir,'AUC_CSdecoder_GLMvsLDA_PCA_',num2str(epoch),'epochs_validatedZ4.jpeg'],'jpeg')
saveas(fh3,[saveDir,'AUC_CSdecoder_GLMvsLDA_PCA_',num2str(epoch),'epochs_validatedZ4.svg'],'svg')
saveas(fh4,[saveDir,'avgAccuracy_CSdecoder_GLMvsLDA_PCA_',num2str(epoch),'epochs_validatedZ4.jpeg'],'jpeg')
saveas(fh4,[saveDir,'avgAccuracy_CSdecoder_GLMvsLDA_PCA_',num2str(epoch),'epochs_validatedZ4.svg'],'svg')
saveas(fh5,[saveDir,'AUC_GLM_CSdecoder_PCA_',num2str(epoch),'epochs_validatedZ4.jpeg'],'jpeg')
saveas(fh5,[saveDir,'AUC_GLM_CSdecoder_PCA_',num2str(epoch),'epochs_validatedZ4.svg'],'svg')
saveas(fh6,[saveDir,'AUC_LDA_CSdecoder_PCA_',num2str(epoch),'epochs_validatedZ4.jpeg'],'jpeg')
saveas(fh6,[saveDir,'AUC_LDA_CSdecoder_PCA_',num2str(epoch),'epochs_validatedZ4.svg'],'svg')

accuracy.GLM = GLMacc;
accuracy.LDA = LDAacc;
accuracy.valGLM = GLMaccShuff;
accuracy.valLDA = LDAaccShuff;
save([dataDir,'Accuracy_CSdecoder_PCA_',num2str(epoch),'epochs_validatedZ4.mat'],'-struct','accuracy');

AUC.GLM = AUCg;
AUC.LDA = AUCl;
AUC.valGLM = AUCgShuff;
AUC.valLDA = AUClShuff;
save([dataDir,'AUC_CSdecoder_PCA_',num2str(epoch),'epochs_validatedZ4.mat'],'-struct','AUC');

close all

save([dataDir,'CSdecoder_VarExplained.mat'],'VarExp')
save([dataDir,'CSdecoder_PCcoeffs.mat'],'PCcoeff')


VarExp = mean(VarExp,3);
sumVarExp = sum(VarExp(1:4,:),1);
semVarExp = std(VarExp')/sqrt(size(VarExp,2));
meanVarExp = mean(VarExp,2);

x = 1:size(meanVarExp,2);
fh0 = figure;
figure(fh0)
plot(sumVarExp,'k-','LineWidth',3)
ylabel('variance explained')
xlabel('Time [sec]')
title('Variance explained by first 4 PCs across epochs')
ylim([50 100])

saveas(fh0,[saveDir,'CSdecoder_PCsVarExplained_acrossTime.jpeg'],'jpeg')
saveas(fh0,[saveDir,'CSdecoder_PCsVarExplained_acrossTime.svg'],'svg')

fh1 = figure;
figure(fh1)
bar(meanVarExp)
hold on
errorbar(meanVarExp,semVarExp)
ylabel('variance explained')
xlabel('PC')

saveas(fh1,[saveDir,'CSdecoder_PCsVarExplained.jpeg'],'jpeg')
saveas(fh1,[saveDir,'CSdecoder_PCsVarExplained.svg'],'svg')

close all