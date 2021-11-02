%% By Neuron - Elastic Net regularized GLM (Logistic) to predict type of US omission based on CS firing (whole CS)

clear all

% get data dir
dataDir = 'C:\Users\LabAdmin\Dropbox\ANALYSIS\Anna-Lena\NAccVTA_April2021\'; % combined Data location
saveDir = 'C:\Users\LabAdmin\Dropbox\ANALYSIS\Anna-Lena\NAccVTA_April2021\images\';

FRdata = load([dataDir,'FR_stage3']);
NClass = load([dataDir,'VTA_values_2SD_stage3.mat']);
NClass = NClass.VTA_values;
Dind = find(cellfun(@(NClass) strcmp(NClass,'d'),NClass(:,5))==1);
CSresponse = load([dataDir,'CSresponse_smoothedDiffScore_2SD_stage3_500ms']);
excInd = strcmp([CSresponse.activity{:,1}],'excited');

% set up data array X (trials x neurons) and y (trials)
type = {'A_nr','B_nr'};

rng default    % Set the seed for reproducibility

epoch = 100;
accuracy = nan(length(Dind),epoch);
accuracyBL = nan(length(Dind),epoch);
relAccuracy = nan(length(Dind),2);
weights = nan(length(Dind),5);
popyTest = [];
popyHat = [];
s = randperm(16);
r = randperm(16);
Eaccuracy = nan(epoch,length(type));
propAccurate = nan(epoch,length(type));
AUC = nan(epoch,1);
indAUC = nan(length(Dind),epoch);
for e = 1:epoch
for d = 1:length(Dind)
    CSData = [];
    USData = [];
    BLdata = [];
    
    for f = 1:length(type)
        if strcmp([type{f}],'A_nr')
            USoff=5;
        elseif strcmp([type{f}],'B_nr')
            USoff=0;
        end
        
        CS = FRdata.(type{f}).VTA.CS_z{Dind(d),1};
        US = FRdata.(type{f}).VTA.US_z{Dind(d),1};
        BL = FRdata.(type{f}).VTA.baseline_z{Dind(d),1};
        
        if strcmp(type{f},'A_nr')
            tempY = ones(size(CS,1),1);
        elseif strcmp(type{f},'B_nr')
            tempY = zeros(size(CS,1),1);
        end
        CS = cat(2,CS(:,e),tempY);
        US = cat(2,US(:,USoff+1:USoff+5),tempY);
        BL = BL(:,e);
        CSData = cat(1,CSData,CS);
        USData = cat(1,USData,US);
        BLdata = cat(1,BLdata,BL);
    end
    % shuffle the row order, so not all trials of the same type follow each
    % other
    BLdata = BLdata(r,:);
    CSData = CSData(s,:);
    yTrain = CSData(:,end);
    CSData = CSData(:,1:end-1);
    USData = USData(s,:);
    yTest = USData(:,end);
    USData = USData(:,1:end-1);
    USData = mean(USData,2);
    
    [B,FitInfo] = lassoglm(USData,yTrain,'binomial','CV',16,'Alpha',.5);
    idxLambdaMinDeviance = FitInfo.IndexMinDeviance;
    B0 = FitInfo.Intercept(idxLambdaMinDeviance);
    coef = [B0; B(:,idxLambdaMinDeviance)];
   
    yhat = glmval(coef,CSData,'logit');
    yHatBinom = (yhat>=0.5);
    
    yhatBL = glmval(coef,BLdata,'logit');
    yHatBinomBL = (yhatBL>=0.5);
    
    [~,~,~,indAUC(d,e)] = perfcurve(yTest,yhat,1);
    
    accuracy(d,e) = mean(yHatBinom == yTest)*100;
    accuracyBL(d,e) = mean(yHatBinomBL == yTest)*100;
    
    cMat(1,1) = length(find(yTest==1 & yHatBinom==1))/length(find(yTest==1)); % 
    cMat(1,2) = length(find(yTest==1 & yHatBinom==0))/length(find(yTest==1));
    cMat(2,1) = length(find(yTest==0 & yHatBinom==1))/length(find(yTest==0));
    cMat(2,2) = length(find(yTest==0 & yHatBinom==0))/length(find(yTest==0));

    popyTest = cat(2,popyTest,yTest);
    popyHat = cat(2,popyHat,yhat);
end

popyTest = mean(popyTest,2);
popyHat = mean(popyHat,2);
[FP,TP,~,AUC(e)] = perfcurve(popyTest,popyHat,1);

close all
end

% identify neurons that perform better or worse than random model
[pop1,~] = find(sum(accuracy(:,3:8)<accuracyBL(:,3:8),2)>3); % population that decodes worse than random
[pop2,~] = find(sum(accuracy(:,3:8)>accuracyBL(:,3:8),2)>3); % population that decodes better than random
pop1 = unique(pop1);
pop2 = unique(pop2);

subpop = struct;
subpop.worseThan = pop1;
subpop.betterThan = pop2;

save([dataDir,'Populations_byNeuron_WorseBetter.mat'],'-struct','subpop');