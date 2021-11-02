function [c,zscore,idxSD]=findOutliersSD(data,column)

% [~,columns] = size(data);
% 
% if columns>1
%     disp('Please only input one column of data for now')
%     c=[];
%     return
% end

smean=mean(data(:,column));
sSD=std(data(:,column));
zscore=nan(length(data(:,column)),1);
SD=zeros(length(data(:,column)),1);
for i=1:length(data(:,column))
    zscore(i)=(data(i,column)-smean)/sSD;
    if zscore(i)>3 ||zscore(i)<-3
        SD(i)=1;
    end
end
idxSD=find(SD==1);
data(idxSD,:)=[];
c=data;