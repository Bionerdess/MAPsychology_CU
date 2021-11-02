%% cross-temporal (US) decoder based on top PCs, compared to LDA
% This protocol runs a PCA for each epoch of CS & period of US omission,
% trains a LDA and GLM on the top PCs of US omission,
% then tests the LDA and GLM on the top PCs of each CS
% period and a GLM the same way.

clear all

% get data dir
dataDir = 'C:\Users\LabAdmin\Dropbox\ANALYSIS\Anna-Lena\NAccVTA_April2021\'; % combined Data location
saveDir = 'C:\Users\LabAdmin\Dropbox\ANALYSIS\Anna-Lena\NAccVTA_April2021\images\';

FRdata = load([dataDir,'FR_stage3']);
NClass = load([dataDir,'VTA_values_2SD_stage3.mat']);
NClass = NClass.VTA_values;
Dind = find(cellfun(@(NClass) strcmp(NClass,'d'),NClass(:,5))==1);

type = {'A_nr','B_nr'};

reps = 1000;
sub = load([dataDir,'Populations_byNeuron_WorseBetterSwitchStay.mat']);
usepops=0; % set to 1 to run multiple decoders for subpopulations, set to 0 to use the whole population
if usepops ==1
    pops = fields(sub);
    
else
    subpop.all = 1:length(Dind);
    pops = {'all'};
    select = 40;
    smpl = randi(length(Dind),reps,select);
end

nBins = 1;
epoch = 100/nBins;
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
perm = nan(reps,16);


VarExp = cell(4,1);
PCAcoeff = cell(4,1);
UScoeffs = struct;
for p = 1:length(pops)
    tic
    if usepops ==1
        select = length(sub.(pops{p}));
        smpl = repmat(sub.(pops{p})',reps,1);
        if select>16
            smpl = randi(select,reps,16);
        end
    end
    for pc = 8
        VarExp{pc} = [];
        PCAcoeff{pc} = [];
        for comb=1:reps
            Y = cat(1,zeros(8,1),ones(8,1)); % trial types (8xA: 0; 8xB: 1)
            s = randperm(16); % shuffle trials, keep the same for all epochs
            Y= Y(s,:); % shuffle Y the same as trial FR data
            perm(comb,:) = s;
            r = randperm(16); % create second shuffle to validate models
            
            USdata = nan(16,select);
            
            score = cell(1,epoch);
            pred = cell(1,epoch);
            
            % set up US PCA, train GLM and LDA
            for du = 1:size(smpl,2)
                USA = FRdata.A_nr.VTA.US_z{Dind(smpl(comb,du)),1};
                USB = FRdata.B_nr.VTA.US_z{Dind(smpl(comb,du)),1};
                
                USdata(:,du) = cat(1,mean(USA(:,6:10),2),mean(USB(:,1:5),2));
            end
            % shuffle US data like Y
            USdata = USdata(s,:);
            
            % run US PCA
            [UScoeff,USscore,~,~,USvarExp] = pca(USdata);
            VarExp{pc} = cat(2,VarExp{pc},USvarExp);
            PCAcoeff{pc} = cat(2,PCAcoeff{pc},UScoeff(:,1:pc));
            
            
            % dot product of top 8 PCs and mean-centered US FR
            UScentered = USdata-mean(USdata,1);
            dotUS = UScentered*UScoeff(:,1:pc);
            % train LDA
            mdl = fitcdiscr(dotUS,Y); % first 8 PCs
            
            % train GLM
            [B,FitInfo] = lassoglm(dotUS,Y,'binomial','CV',8,'Alpha',.5); % train LASSO GLM on top 8 PCA scores
            idxLambdaMinDeviance = FitInfo.IndexMinDeviance;
            B0 = FitInfo.Intercept(idxLambdaMinDeviance);
            coef = [B0; B(:,idxLambdaMinDeviance)];
            
            for e=1:epoch
                neuronData = [];
                BLdata = [];
                for f = 1:length(type)
                    temp = nan(8,select);
                    tempBL = nan(8,select);
                    for d = 1:select
                        CS = FRdata.(type{f}).VTA.CS_z{Dind(smpl(comb,d)),1};
                        BL = FRdata.(type{f}).VTA.baseline_z{Dind(smpl(comb,d)),1};
                        temp(:,d) = mean(CS(:,e*nBins-(nBins-1):e*nBins),2);
                        tempBL(:,d) = mean(BL(:,e*nBins-(nBins-1):e*nBins),2);
                    end
                    
                    neuronData = cat(1,neuronData,temp);
                    BLdata = cat(1,BLdata,tempBL);
                end
                % shuffle trials
                neuronData = neuronData(s,:);
                centeredDat = neuronData-mean(neuronData,1);
                BLdata = BLdata(s,:);
                centeredBL = BLdata-mean(BLdata,1);
                
                dotDat = centeredDat*UScoeff(:,1:pc);
                dotShuff = centeredBL*UScoeff(:,1:pc);
                
                % test LDA on CS epoch
                pred{e} = predict(mdl,dotDat);
                c = confusionmat(Y,pred{e});
                [~,~,~,AUCl(comb,e)] = perfcurve(Y,pred{e},1);
                
                LDAacc(comb,e,1) = c(1,1)/(c(1,1)+c(1,2)); % percent accurate out of A trials
                LDAacc(comb,e,2) = c(2,2)/(c(2,1)+c(2,2));
                
                % test validation LDA on shuffled data
                predShuff = predict(mdl,dotShuff);
                cShuff = confusionmat(Y,predShuff);
                [~,~,~,AUClShuff(comb,e)] = perfcurve(Y,predShuff,1);
                
                LDAaccShuff(comb,e,1) = cShuff(1,1)/(cShuff(1,1)+cShuff(1,2)); % percent accurate out of A trials
                LDAaccShuff(comb,e,2) = cShuff(2,2)/(cShuff(2,1)+cShuff(2,2));
                
                
                % test GLM on CS epoch
                yhat = glmval(coef,dotDat,'logit');
                yHatBinom = double(yhat>=0.5);
                
                [~,~,~,AUCg(comb,e)] = perfcurve(Y,yhat,1);
                
                [mat,gorder] = confusionmat(Y,yHatBinom);
                
                GLMacc(comb,e,1) = mat(1,1)/(mat(1,1)+mat(1,2));
                GLMacc(comb,e,2) = mat(2,2)/(mat(2,1)+mat(2,2));
                
                % test validation GLM on shuffled Data
                yhatShuff = glmval(coef,dotShuff,'logit');
                yHatBinomShuff = double(yhatShuff>=0.5);
                
                [~,~,~,AUCgShuff(comb,e)] = perfcurve(Y,yhatShuff,1);
                
                matShuff = confusionmat(Y,yHatBinomShuff);
                
                GLMaccShuff(comb,e,1) = matShuff(1,1)/(matShuff(1,1)+matShuff(1,2));
                GLMaccShuff(comb,e,2) = matShuff(2,2)/(matShuff(2,1)+matShuff(2,2));
                
            end
            
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
        title([pops{p},' - n= ',num2str(select)])
        xlabel('Time (100ms bins)')
        ylabel('% accurate')
        
        fh2 = figure;
        figure(fh2)
        plot(mean(mLDAacc,2),'LineWidth',2)
        hold on
        plot(mean(mLDAaccShuff,2),'LineWidth',2)
        hold on
        plot([0 epoch],[.5 .5],'k--')
        ylim([0 1])
        legend('LDAacc','validation LDAacc')
        title([pops{p},' - n= ',num2str(select)])
        xlabel('Time (100ms bins)')
        ylabel('% accurate')
        
        fh3 = figure;
        figure(fh3)
        plot(mean(AUCg,1),'LineWidth',3)
        hold on
        plot(mean(AUCl,1),'Linewidth',3)
        hold on
        plot([0 epoch],[.5 .5],'k--')
        ylim([0 1])
        legend('GLM','LDA')
        title([pops{p},' - n= ',num2str(select)])
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
        title([pops{p},' - n= ',num2str(select)])
        xlabel('Time (100ms bins)')
        ylabel('avg accuracy')
        
        fh5 = figure;
        figure(fh5)
        plot(mean(AUCg,1),'LineWidth',3)
        hold on
        plot(mean(AUCgShuff,1),'Linewidth',3)
        hold on
        plot([0 epoch],[.5 .5],'k--')
        ylim([0 1])
        legend('GLM','validation GLM')
        title([pops{p},' - n= ',num2str(select)])
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
        title([pops{p},' - n= ',num2str(select)])
        xlabel('Time (100ms bins)')
        ylabel('AUC')
        
        saveas(fh1,[saveDir,'avgAccuracyGLM_USdecoder_PC',num2str(pc),'_',pops{p},'_',num2str(epoch),'epochs_validatedZ4.jpeg'],'jpeg')
        saveas(fh1,[saveDir,'avgAccuracyGLM_USdecoder_PC',num2str(pc),'_',pops{p},'_',num2str(epoch),'epochs_validatedZ4.svg'],'svg')
        saveas(fh2,[saveDir,'avgAccuracyLDA_USdecoder_PC',num2str(pc),'_',pops{p},'_',num2str(epoch),'epochs_validatedZ4.jpeg'],'jpeg')
        saveas(fh2,[saveDir,'avgAccuracyLDA_USdecoder_PC',num2str(pc),'_',pops{p},'_',num2str(epoch),'epochs_validatedZ4.svg'],'svg')
        saveas(fh3,[saveDir,'AUC_USdecoder_PC',num2str(pc),'_',pops{p},'_GLMvsLDA_',num2str(epoch),'epochs_validatedZ4.jpeg'],'jpeg')
        saveas(fh3,[saveDir,'AUC_USdecoder_PC',num2str(pc),'_',pops{p},'_GLMvsLDA_',num2str(epoch),'epochs_validatedZ4.svg'],'svg')
        saveas(fh4,[saveDir,'avgAccuracy_USdecoder_PC',num2str(pc),'_',pops{p},'_GLMvsLDA_',num2str(epoch),'epochs_validatedZ4.jpeg'],'jpeg')
        saveas(fh4,[saveDir,'avgAccuracy_USdecoder_PC',num2str(pc),'_',pops{p},'_GLMvsLDA_',num2str(epoch),'epochs_validatedZ4.svg'],'svg')
        saveas(fh5,[saveDir,'AUC_GLM_USdecoder_PC',num2str(pc),'_',pops{p},'_',num2str(epoch),'epochs_validatedZ4.jpeg'],'jpeg')
        saveas(fh5,[saveDir,'AUC_GLM_USdecoder_PC',num2str(pc),'_',pops{p},'_',num2str(epoch),'epochs_validatedZ4.svg'],'svg')
        saveas(fh6,[saveDir,'AUC_LDA_USdecoder_PC',num2str(pc),'_',pops{p},'_',num2str(epoch),'epochs_validatedZ4.jpeg'],'jpeg')
        saveas(fh6,[saveDir,'AUC_LDA_USdecoder_PC',num2str(pc),'_',pops{p},'_',num2str(epoch),'epochs_validatedZ4.svg'],'svg')
        
        close all
        
        accuracy.GLM = GLMacc;
        accuracy.LDA = LDAacc;
        accuracy.valGLM = GLMaccShuff;
        accuracy.valLDA = LDAaccShuff;
        save([dataDir,'Accuracy_USdecoder_PC',num2str(pc),'_',pops{p},'_',num2str(epoch),'epochs_validatedZ4.mat'],'-struct','accuracy');
        
        AUC.GLM = AUCg;
        AUC.LDA = AUCl;
        AUC.valGLM = AUCgShuff;
        AUC.valLDA = AUClShuff;
        save([dataDir,'AUC_USdecoder_PC',num2str(pc),'_',pops{p},'_',num2str(epoch),'epochs_validatedZ4.mat'],'-struct','AUC');
        toc
    end
end

UScoeffs.PC4 = PCAcoeff{4};
UScoeffs.smpl = smpl';

save([dataDir,'USdecoder_coeff'],'-struct','UScoeffs')

VarExp = VarExp{pc};
meanVarExp = mean(VarExp,2);
semVarExp = std(VarExp')/sqrt(size(VarExp,2));
fh0 = figure;
figure(fh0)
bar(meanVarExp)
hold on
errorbar(meanVarExp,semVarExp)
xlim([0 16])
ylabel('variance explained')
xlabel('PC')

saveas(fh0,[saveDir,'USdecoder_PCsVarExplained.jpeg'],'jpeg')
saveas(fh0,[saveDir,'USdecoder_PCsVarExplained.svg'],'svg')
save([dataDir,'USdecoder_PCs_VarExplained.mat'],'VarExp')