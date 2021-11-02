function  varargout = ClusterAnalysis(x,k,rep)

%ClusterAnalysis - exploratory cluster analysis using kmeans function
%
%   IDX = ClusterAnalysis(x,k,rep) runs kmeans function k 
%   times for 1:k number of clusters on N-by-P matrix X. 
%   Rows (N) of dat correspond to points, while Columns (P) correspond to 
%   variables. Each iteration will run rep number of repetitions. IDX is a 
%   M-by-N matrix with M providing the cluster affiliation for each point.
%   This holds true for each Nth column corresponding to N number of
%   clusters utilized to calculate kmeans.
%   
%   [IDX,ED] = ClusterAnalysis(x,k,rep) returns the mean squared
%   Euclidean distance for each iteration of k to provide information about
%   the optimal number of clusters.
%
%   [IDX,ED,FIG] = ClusterAnalysis(x,k,rep) returns a bar graph of the
%   mean squared Euclidean distance in regards to number of clusters to
%   provide information about the optimal number of clusters.
%
%   Note: This cluster analysis uses the Euclidean distance as distance
%   measure. For further documentation on the kmeans function see KMEANS
%
%
%   Anna-Lena Schlenner
%   al.schlenner@gmail.com
%
%   April 2021

[N,P] = size(x);
sc = nan(N,k);
hist = nan(1,k);

for n=1:k
    [ID,~,sumD] = kmeans(x,n,'Replicates',rep);
    sc(:,n) = ID;
    hist(n) = mean(sumD);
end

fh1=figure;
figure(fh1)
bar(hist)
xlabel('nCluster')
ylabel('Mean Squared Euclidean Distance')
title('Mean Squared Euclidean Distance across number of clusters')

varargout{1} = sc;
varargout{2} = hist;
varargout{3} = fh1;