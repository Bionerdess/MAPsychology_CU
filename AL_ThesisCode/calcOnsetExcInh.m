%% define excitation & inhibition onset

clear all

% data Directory

dataDir = 'C:\Users\LabAdmin\Dropbox\ANALYSIS\Anna-Lena\NAccVTA_April2021\'; % combined Data location
saveDir = 'C:\Users\LabAdmin\Dropbox\ANALYSIS\Anna-Lena\NAccVTA_April2021\images\';

nSig = 3;
window = 97;

statsCS = load([dataDir,'ExcInhPeriods_AB_CS.mat']);
statsCS = statsCS.CSstat;

res = cell(2,size(statsCS,2));
Binh = zeros(length(statsCS),100);
Aexc = zeros(length(statsCS),100);
Bexc = zeros(length(statsCS),100);
Ainh = zeros(length(statsCS),100);

s = randperm(100);

for c = 1:size(statsCS,2)
    Asum = nan(1,window);
    Bsum = nan(1,window);
    
    for i = 1:window
        Asum(1,i) = sum(statsCS{1,c}(i:i+nSig-1));
        Bsum(1,i) = sum(statsCS{2,c}(i:i+nSig-1));
    end
    res{1,c}{1} = find(Asum ==nSig);
    res{1,c}{2} = find(Asum ==-nSig);
    res{2,c}{1} = find(Bsum ==nSig);
    res{2,c}{2} = find(Bsum ==-nSig);
    
    if ~isempty(res{2,c}{2})
        Binh(c,res{2,c}{2}(1)) = 1;
    end
    if ~isempty(res{1,c}{1})
        Aexc(c,res{1,c}{1}(1)) = 1;
    end
    if ~isempty(res{2,c}{1})
        Bexc(c,res{2,c}{1}(1)) = 1;
    end
    if ~isempty(res{1,c}{2})
        Ainh(c,res{1,c}{2}(1)) = 1;
    end
end

Exc(1,:) = sum(Aexc,1);
Exc(2,:) = sum(Bexc,1);
Inh(1,:) = sum(Ainh,1);
Inh(2,:) = sum(Binh,1);

shuffE(1,:) = sum(Aexc(:,s),1);
shuffE(2,:) = sum(Bexc(:,s),1);
shuffI(1,:) = sum(Ainh(:,s),1);
shuffI(2,:) = sum(Binh(:,s),1);

[statsE(1,:),cumSumE(1,:),ShuffCumSumE(1,:)] = calcCumSum(Exc(1,:),shuffE(1,:),'alpha',0.001);
[statsE(2,:),cumSumE(2,:),ShuffCumSumE(2,:)] = calcCumSum(Exc(2,:),shuffE(2,:),'alpha',0.001);
[statsI(1,:),cumSumI(1,:),ShuffCumSumI(1,:)] = calcCumSum(Inh(1,:),shuffI(1,:),'alpha',0.001);
[statsI(2,:),cumSumI(2,:),ShuffCumSumI(2,:)] = calcCumSum(Inh(2,:),shuffI(2,:),'alpha',0.001);

ExcAU = AUROC(Exc(1,:),Exc(2,:));
InhAU = AUROC(Inh(1,:),Inh(2,:));

excAU = AUROC(cumSumE(1,:),cumSumE(2,:));
inhAU = AUROC(cumSumI(1,:),cumSumE(2,:));

[hE,pE,statE] = kstest2(cumSumE(1,:),cumSumE(2,:));
[hI,pI,statI] = kstest2(cumSumI(1,:),cumSumI(2,:));

fh1 = figure;
figure(fh1)
plot(cumSumE(1,:),'r')
hold on
plot(cumSumE(2,:),'k')
ylabel('n Neurons')
xlabel('Time 100ms bins')
title('Excitation onset')
legend('A','B','location','southeast')

fh3 = figure;
figure(fh3)
plot(cumSumI(1,:),'r')
hold on
plot(cumSumI(2,:),'k')
ylabel('n Neurons')
xlabel('Time 100ms bins')
title('Inhibition onset')
legend('A','B','location','southeast')


saveas(fh1,[saveDir,'ExcInh_Onset_ExcAvsB.jpeg'],'jpeg')
saveas(fh3,[saveDir,'ExcInh_Onset_InhAvsB.jpeg'],'jpeg')
