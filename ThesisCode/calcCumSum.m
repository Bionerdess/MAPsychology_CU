function [stats,cumSumSmpl,cumSumBtstrp] = calcCumSum(smpl,btstrpSmpl,varargin)

%%% calcCumSum takes an M x N sample matrix and compares the cumulative 
%%% sum of the sample to a bootstrapped "baseline" sample. 
%%% Both matrixes have to be shaped as observations x time with the same
%%% number of columns.
%%%
%%% [stats] = calcCumSum(...) returns an array identifying whether the
%%% change in cumulative sum of the sample is significantly smaller (-1), 
%%% larger (1) or not different (0) from the "baseline" sample.
%%% 
%%% [...,cumSumSmpl,cumSumBtstrp] = calcCumSum(...) returns the cumulative
%%% sum of the sample and bootstrapped "baseline" sample respectively.
%%%
%%% [...] = calcCumSum(...,'nBoot',nBoot) specifies the number of
%%% bootstrapped samples (Default: 200)
%%%
%%% [...] = calcCumSum(...,'alpha',alpha) uses the specified alpha value to
%%% calculate the bootstrapped confidence interval. (Default: alpha = 0.05)
%%% 
%%% [...] = calcCumSum(...,'tail',tail) specifies one-sided (1) or
%%% two-sided (2) interval bounds (Default: tail = 2)
%%%
%%%
%%% Anna-Lena Schlenner, August 2021
%%% al.schlenner@gmail.com

if size(smpl,2)~= size(btstrpSmpl,2)
    error(message('signal:calcCumSum:MatrixDimensionMismatch', num2str(size(smpl,2)),num2str(size(btstrpSmpl,2))));
end

if ~isempty(varargin)
    [alpha,tail,nBoot] = getargs(varargin);
else
    alpha = 0.05;
    tail = 2;
    nBoot = 200;
end

% calc cumsum and plot
avgB = mean(btstrpSmpl,1);
avgS = mean(smpl,1);

CI = nan(2,size(btstrpSmpl,2));
stats = zeros(1,size(btstrpSmpl,2));
% bootstrap BL data & calc cumsum (& shuffle?)
[bootstat,~] = bootstrp(nBoot,@cumsum,avgB');

% calc CI for each sample in bootstat

cumSumBtstrp = mean(bootstat,1);
% df = size(bootstat,1)-1;

% sort each column small to large, (alpha/tail)* size(~,1)
% count that many values from bottom and top to
% identify lower and upper bound of CI
nVal = floor((alpha/tail)*size(bootstat,1));
for z=1:size(bootstat,2)
    bootstat(:,z) = sort(bootstat(:,z));
    CI(1,z) = bootstat(1+nVal,z);
    CI(2,z) = bootstat(size(bootstat,1)-nVal,z);
end

cumSumSmpl = cumsum(avgS);
diffS = diff(cumSumSmpl); % same as avgCS(2:end)
diffCI = diff(CI,1,2);

for q=1:length(diffS)
    
    if diffS(q)<diffCI(1,q) % if less than lower CI bound
        stats(1,q) = -1;
    elseif diffS(q) > diffCI(2,q) % if greater than upper CI bound
        stats(1,q) = 1;
    end
end
                        
%---------------------------------------                        
                        
function [alpha,tail,nBoot] = getargs(args)
p = inputParser;
for a=1:2:size(args,2)
    addParameter(p,args{a},[]);
end

parse(p,args{:});
alpha = 0.05;
tail=2;
nBoot = 200;
for b=1:length(p.Parameters)
    if strcmp(p.Parameters{b},'alpha')
        alpha = p.Results.(p.Parameters{b});
    end
    if strcmp(p.Parameters{b},'nBoot')
        nBoot = p.Results.(p.Parameters{b});
    end
    if strcmp(p.Parameters{b},'tail')
        tail = p.Results.(p.Parameters{b});
        if ischar(tail)
            error('Error. \nInput must be 1 or 2, not a %s.',class(tail))
        elseif tail~=1 && tail~=2
            error('Error. \nInput must be 1 or 2, not  %s.',num2str(tail))
        end
    end
end