function [tSp,tSp_bin,good_neuron_ids,peakCh] = gettSp_bin(dir)
%%% identify spike times of neurons (based on phyed clusters)
%%% determine neuron IDs, peak channel etc.

cd(dir)

% read Phy Results nad unfold neurons
addpath(genpath('C:\Users\LabAdmin\Dropbox\KiloSort')) % path to kilosort folder

% Import CSV cluster data

filename = 'cluster_info.tsv';

delimiter = '\t';
startRow = 2;
formatSpec = '%f%f%f%f%s%s%s%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);

cluster_id = dataArray{:, 1};
cluster_type = dataArray{:, 6};
channel = dataArray{:,3};

% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

% initalize unfolded neuron table

 SPIKE.times = readNPY('spike_times.npy');
 SPIKE.clusters = readNPY('spike_clusters.npy');

k=1; 
for i=1:length(cluster_type) %get IDS of good neurons
    if(strcmp(cluster_type(i),'good'))
        good_neuron_ids(k) =  cluster_id(i);
        peakCh(k) = channel(i);
        k=k+1;
    end
end
count=0;
if exist('good_neuron_ids','var')
    for i=1:length(good_neuron_ids)
        for r=1:length(SPIKE.clusters)
            
            if(good_neuron_ids(i)==SPIKE.clusters(r))
                
                EXPANDED_Neurons(i,SPIKE.times(r))=1;
                count=count+1;
            end
            
        end
    end
    
    [nCells,~] = size(EXPANDED_Neurons);
    for n=1:nCells
        tSp{n} = find(EXPANDED_Neurons(n,:)==1);
    end
else
    tSp=[];
    good_neuron_ids = [];
    peakCh = [];
    disp('No good clusters!')
end
tSp_bin = EXPANDED_Neurons;