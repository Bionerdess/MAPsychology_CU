function [MainA,MainB,Interaction,SimpleAatB,SimpleBatA] = mik_anova2(F1,F2,S,Y)
% called with mystats.m

stats = cell(4,5);

F1_lvls = unique(F1);
F2_lvls = unique(F2);
Subjs = unique(S);

a = length(F1_lvls); % # of levels in factor 1
b = length(F2_lvls); % # of levels in factor 2
n = length(Subjs); % # of subjects

INDS = cell(a,b,n); % this will hold arrays of indices
CELLS = cell(a,b,n); % this will hold the data for each subject X condition
MEANS = zeros(a,b,n); % this will hold the means for each subj X condition

% Calculate means for each subject X condition.
% Keep data in CELLS, because in future we may want to allow options for
% how to compute the means (e.g. leaving out outliers > 3stdev, etc...).
for i=1:a % F1
    for j=1:b % F2
        for k=1:n % Subjs
            INDS{i,j,k} = find(F1==F1_lvls(i) & F2==F2_lvls(j) & S==Subjs(k));
            CELLS{i,j,k} = Y(INDS{i,j,k});
            MEANS(i,j,k) = mean(CELLS{i,j,k});
        end
    end
end

% make tables (see table 18.1, p. 402)
AB = reshape(sum(MEANS,3),a,b); % across subjects
AS = reshape(sum(MEANS,2),a,n); % across factor 2
BS = reshape(sum(MEANS,1),b,n); % across factor 1

A = sum(AB,2); % sum across columns, so result is ax1 column vector
B = sum(AB,1); % sum across rows, so result is 1xb row vector
S = sum(AS,1); % sum across columns, so result is 1xs row vector
T = sum(sum(A)); % could sum either A or B or S, choice is arbitrary

% degrees of freedom
dfA = a-1;
dfB = b-1;
dfAB = (a-1)*(b-1);
dfS = n-1;
dfAS = (a-1)*(n-1);
dfBS = (b-1)*(n-1);
dfABS = (a-1)*(b-1)*(n-1);

% bracket terms (expected value)
expA = sum(A.^2)./(b*n);
expB = sum(B.^2)./(a*n);
expAB = sum(sum(AB.^2))./n;
expS = sum(S.^2)./(a*b);
expAS = sum(sum(AS.^2))./b;
expBS = sum(sum(BS.^2))./a;
expY = sum(Y.^2);
expT = T^2 / (a*b*n); % (otherwise named CF in Hyperstats)

% sums of squares
ssA = expA - expT;
ssB = expB - expT;
ssAB = expAB - expA - expB + expT;
ssS = expS - expT;
ssAS = expAS - expA - expS + expT;
ssBS = expBS - expB - expS + expT;
ssABS = expY - expAB - expAS - expBS + expA + expB + expS - expT;
ssTot = expY - expT;

% mean squares
msA = ssA / dfA;
msB = ssB / dfB;
msAB = ssAB / dfAB;
msS = ssS / dfS;
msAS = ssAS / dfAS;
msBS = ssBS / dfBS;
msABS = ssABS / dfABS;

% f statistic
fA = msA / msAS;
fB = msB / msBS;
fAB = msAB / msABS;

% p values
pA = 1-fcdf(fA,dfA,dfAS);
pB = 1-fcdf(fB,dfB,dfBS);
pAB = 1-fcdf(fAB,dfAB,dfABS);
disp(' ');disp('Two-way ANOVA:');
disp(['Main effect A: F(' num2str(dfA) ',' num2str(dfAS) ') = ' num2str(fA) ', p = ' num2str(pA) '.']);
disp(['Main effect B: F(' num2str(dfB) ',' num2str(dfBS) ') = ' num2str(fB) ', p = ' num2str(pB) '.']);
disp(['2-way interaction A-B: F(' num2str(dfAB) ',' num2str(dfABS) ') = ' num2str(fAB) ', p = ' num2str(pAB) '.']);

MainA = table(pA,fA,dfA,dfAS,'VariableNames',{'pA','FA','df','dfError'});
MainB = table(pB,fB,dfB,dfBS,'VariableNames',{'pB','FB','df','dfError'});
Interaction = table(pAB,fAB,dfAB,dfABS,'VariableNames',{'pAB','fAB','df','dfError'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simple main effects on 2-ways interactions
expA_atBj=[];expB_atAi=[];
for i=1:a
    expB_atAi(i) = sum(AB(i,:).^2)/n;
end
for j=1:b
    expA_atBj(j) = sum(AB(:,j).^2)/n;
end

% sums of squares
ssA_atBj = expA_atBj - (B.^2)/(a*n);
ssB_atAi = expB_atAi - (A'.^2)/(b*n);

% mean squares
msA_atBj = ssA_atBj / dfA;
msB_atAi = ssB_atAi / dfB;

% f statistic
fA_atBj = msA_atBj ./ msAS;
fB_atAi = msB_atAi ./ msBS;

AatB = nan(b,4);
BatA = nan(b,4);
disp(' ');disp('Simple main effects across 0 factor:');
for j=1:b
    disp(['Simple effect of A at levelB' num2str(j) ': F(' num2str(dfA) ',' num2str(dfAS) ') = ' num2str(fA_atBj(j)) ', p = ' num2str(1-fcdf(fA_atBj(j),dfA,dfAS)) '.']);
    AatB(j,:) = [1-fcdf(fA_atBj(j),dfA,dfAS),fA_atBj(j),dfA,dfAS];
end
for i=1:a
    disp(['Simple effect of B at levelA' num2str(i) ': F(' num2str(dfB) ',' num2str(dfBS) ') = ' num2str(fB_atAi(i)) ', p = ' num2str(1-fcdf(fB_atAi(i),dfB,dfBS)) '.']);
    BatA(i,:) = [1-fcdf(fB_atAi(i),dfB,dfBS),fB_atAi(i),dfB,dfBS];
end

SimpleAatB = table(AatB(:,1),AatB(:,2),AatB(:,3),AatB(:,4),'VariableNames',{'pAatB','FAatB','df','dfError'});
SimpleBatA = table(BatA(:,1),BatA(:,2),BatA(:,3),BatA(:,4),'VariableNames',{'pBatA','FBatA','df','dfError'});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MB=B/(a*n);
MA=A/(b*n);
critqA=stu_dist(dfAS,dfA);
critqB=stu_dist(dfBS,dfB);

qA=[];qB=[];
disp(' ');disp('TUKEY''s HSD Pairwise comparisons between levels of factorA:');
for i=1:a
    for co=i+1:a

        if dfA==1
            qA(co,i)=abs(MA(i)-MA(co))/sqrt(2*msAS/(b*n));h=floor(qA(co,i)/critqA);
            disp('With only 2 means, there are no multiple comparisons. So Tukey is basically a ttest, i.e. sqrt(F value).');
            if h>=1;disp(['SIGNIFICANT difference between levelA' num2str(i) ' and levelA' num2str(co) ': q = ' num2str(qA(co,i)) ', p=' num2str(1-fcdf(qA(co,i)^2,dfA,dfAS)) '.']);
            else disp(['non-sign. difference between levelA' num2str(i) ' and levelA' num2str(co) ': q = ' num2str(qA(co,i)) ', p=' num2str(1-fcdf(qA(co,i)^2,dfA,dfAS)) '.']);end
        else
            qA(co,i)=abs(MA(i)-MA(co))/sqrt(msAS/(b*n));h=floor(qA(co,i)/critqA);
            if h>=1;disp(['SIGNIFICANT difference between levelA' num2str(i) ' and levelA' num2str(co) ': q = ' num2str(qA(co,i)) ', p<0.05.']);
            else disp(['non-sign. difference between levelA' num2str(i) ' and levelA' num2str(co) ': q = ' num2str(qA(co,i)) ', p>0.05.']);end
        end
        
    end
end
disp(' ');disp('TUKEY''s HSD Pairwise comparisons between levels of factorB:');
for i=1:b
    for co=i+1:b
        
        if dfB==1
            qB(co,i)=abs(MB(i)-MB(co))/sqrt(2*msBS/(a*n));h=floor(qB(co,i)/critqB);
            disp('With only 2 means, there are no multiple comparisons. So Tukey is basically a ttest, i.e. sqrt(F value).');
            if h>=1;disp(['SIGNIFICANT difference between levelB' num2str(i) ' and levelB' num2str(co) ': q = ' num2str(qB(co,i)) ', p=' num2str(1-fcdf(qB(co,i)^2,dfB,dfBS)) '.']);
            else disp(['non-sign. difference between levelB' num2str(i) ' and levelB' num2str(co) ': q = ' num2str(qB(co,i)) ', p=' num2str(1-fcdf(qB(co,i)^2,dfB,dfBS)) '.']);end
        else
            qB(co,i)=abs(MB(i)-MB(co))/sqrt(msBS/(a*n));h=floor(qB(co,i)/critqB);
            if h>=1;disp(['SIGNIFICANT difference between levelB' num2str(i) ' and levelB' num2str(co) ': q = ' num2str(qB(co,i)) ', p<0.05.']);
            else disp(['non-sign. difference between levelB' num2str(i) ' and levelB' num2str(co) ': q = ' num2str(qB(co,i)) ', p>0.05.']);end
        end
        
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('  ');disp('Tukey''s HSD pairwise comparisons across 0 factor:');disp(' ');

MAB=AB/n;
qA_at_Bi=[];qB_at_Ai=[];

disp(' ');disp('TUKEY''s HSD Pairwise comparisons between levels of factorA at each level of B:');disp(' ');
for j=1:b
    disp(['LevelB' num2str(j) ]);
    for i=1:a
        for co=i+1:a

            if dfA==1
                qA_at_Bi(co,i)=abs(MAB(i,j)-MAB(co,j))/sqrt(2*msAS/n);h=floor(qA_at_Bi(co,i)/critqA);
                if h>=1;disp(['SIGNIFICANT difference between levelA' num2str(i) ' and levelA' num2str(co) ': q = ' num2str(qA_at_Bi(co,i)) ', p=' num2str(1-fcdf(qA_at_Bi(co,i)^2,dfA,dfAS)) '.']);
                else disp(['non-sign. difference between levelA' num2str(i) ' and levelA' num2str(co) ': q = ' num2str(qA_at_Bi(co,i)) ', p=' num2str(1-fcdf(qA_at_Bi(co,i)^2,dfA,dfAS)) '.']);end
            else
                qA_at_Bi(co,i)=abs(MAB(i,j)-MAB(co,j))/sqrt(msAS/n);h=floor(qA_at_Bi(co,i)/critqA);
                if h>=1;disp(['SIGNIFICANT difference between levelA' num2str(i) ' and levelA' num2str(co) ': q = ' num2str(qA_at_Bi(co,i)) ', p<0.05.']);
                else disp(['non-sign. difference between levelA' num2str(i) ' and levelA' num2str(co) ': q = ' num2str(qA_at_Bi(co,i)) ', p>0.05.']);end
            end
            
        end
    end
end
disp(' ');disp('TUKEY''s HSD Pairwise comparisons between levels of factorB at each level of A:');disp(' ');
for j=1:a
    disp(['LevelA' num2str(j) ]);
    for i=1:b
        for co=i+1:b
            
            if dfB==1
                qB_at_Ai(co,i)=abs(MAB(j,i)-MAB(j,co))/sqrt(2*msBS/n);h=floor(qB_at_Ai(co,i)/critqB);
                if h>=1;disp(['SIGNIFICANT difference between levelB' num2str(i) ' and levelB' num2str(co) ': q = ' num2str(qB_at_Ai(co,i)) ', p=' num2str(1-fcdf(qB_at_Ai(co,i)^2,dfB,dfBS)) '.']);
                else disp(['non-sign. difference between levelB' num2str(i) ' and levelB' num2str(co) ': q = ' num2str(qB_at_Ai(co,i)) ', p=' num2str(1-fcdf(qB_at_Ai(co,i)^2,dfB,dfBS)) '.']);end
            else
                qB_at_Ai(co,i)=abs(MAB(j,i)-MAB(j,co))/sqrt(msBS/n);h=floor(qB_at_Ai(co,i)/critqB);
                if h>=1;disp(['SIGNIFICANT difference between levelB' num2str(i) ' and levelB' num2str(co) ': q = ' num2str(qB_at_Ai(co,i)) ', p<0.05.']);
                else disp(['non-sign. difference between levelB' num2str(i) ' and levelB' num2str(co) ': q = ' num2str(qB_at_Ai(co,i)) ', p>0.05.']);end
            end
            
        end
    end
end

% disp(' ');disp(' ');
% disp(['To get the exact p values for factorA, refer to: onlinestatbook.com/calculators/tukeycdf.html with Nb_means=' num2str(dfA) ' and df=' num2str(dfAS) '.']);
% disp(['To get the exact p values for factorB, refer to: onlinestatbook.com/calculators/tukeycdf.html with Nb_means=' num2str(dfB) ' and df=' num2str(dfBS) '.']);
% disp(' ');disp(' ');


