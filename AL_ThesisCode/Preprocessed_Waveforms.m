% *************************************************************************
% This script will produce average waveforms, using spike times and clusers
% generated by Kilo/Phy, and will apply kilo stlye preproccessing to raw
% headstage data from the WhiteMatter system

% stephen.gv@gmail.com
% 7/20/18
%
% Iordanova Lab Edits by Anna-Lena Schlenner
% al.schlenner@gmail.com
% 4/21/2020
% *************************************************************************


% EDIT: 
% bandpass only peakCh of good clusters
% clear hs, Median_hs... as soon as not needed anymore!!!
% auto loop through all sessions
% save AvgSpike & SEMSpike variables for all sessions
% 
% 


% **********************PARAMETERS*****************************************
clear all
nchans = 64; % number of channels in headstage data

[b1, a1] = butter(3, .024, 'high'); % Highpass filter creation (.024 = 300 hz high pass, expressed in proportion of nyquist

animalData = load('C:\Users\LabAdmin\Dropbox\ANALYSIS\Anna-Lena\animalDataNAccVTA'); %Load the Data Def
path = animalData.path;
nSessions = length(path);
type = 'NAccVTA_fear';

clear animalData

tic
% **********************LOAD DATA***************************************************************

for j = 67 % :nSessions % 3/27/21 1:40pm
    
    dataDir = path{j};
    cd (dataDir)
    disp(['Processing data of ',dataDir])
    
    saveDir = [dataDir,'\Analysis_April2021'];
    
    % Load cluster and spike data
    %
    
    rez=load('rez.mat');
    rez=rez.rez;
    rez=rez.WrotInv;
    'Whitening Matrix loaded' %#ok<NOPTS>
    toc
    
    %
    
    filehs = dir([dataDir '\H_NOHEADER*.bin']);
    filehs = filehs.name;
    hs = []; %extracted continous data
    ns = int64([]); % number of samples in each headstage file
    fid = fopen([dataDir '\' filehs],'r');
    hs = fread(fid,[nchans,inf],'int16=>int16');
    ns=length(hs);
    'Raw data loaded' %#ok<NOPTS>
    toc
    
    
% **********************APPLY PREPROCESSING***************************************************************

    tic
    Median_hs = zeros(nchans,ns);
    for b=1:ns
        Median_hs(:,b) = hs(:,b) - median(hs(:,b));
    end
    'Common Average Reference Finished' %#ok<NOPTS>
    toc
    clear hs
    
     % Load peakCh
    cellInfo = load([dataDir,'\Analysis_April2021\','CellInfo_',type,'.mat']);
    peakCh = cellInfo.peakCh;
    peakCh=[peakCh{:}]+1; % because phy output 0-63 while chanMap 1-64
    uniqueID = cellInfo.uniqueID;
    
    clear cellInfo
    
        % White Matter EIB Channel Map
    peakCh_WMMap = nan(length(peakCh),1);
    WM_ChanMap = load('C:\Users\LabAdmin\Dropbox\KiloSort\ChannelMap.mat');
    for c=1:length(peakCh)
        peakCh_WMMap(c) = WM_ChanMap.chanMap(peakCh(c)); %fix!!!!
    end

    
    Filt_hs = zeros(length(peakCh_WMMap),ns);
    for a=1:length(peakCh_WMMap) %create bandpassed contiin signal
        
        % Filt_hs(i,:)= filter(butt,double(hs(i,:)));
        Filt_hs(a,:)=filtfilt(b1,a1,double(Median_hs(peakCh_WMMap(a),:)));
        
        ['Channel ' int2str(peakCh_WMMap(a)) ' Highpassed']
        toc
    end
    'Highpass Finished' %#ok<NOPTS>
    toc
    
    clear Median_hs
    plot(Filt_hs(1,1:1000000))
    
    %
    
%     [nChan,~]= size(Filt_hs);
%     Median_hs=zeros(length(peakCh_WMMap),ns);
%     if nChan>1
%         for b=1:ns
%             Median_hs(:,b) = Filt_hs(:,b) - median(Filt_hs(:,b));
%         end
%     else
%         Median_hs(:) = Filt_hs(:);
%     end
    

    %clear Filt_hs
%     plot(Median_hs(1,1:1000000))
    
    tic
    for c=1:length(peakCh_WMMap)
        White_hs(c,:)=Filt_hs(c,:)'*rez(peakCh_WMMap(c));
        %     White_hs=White_hs';
        'Data Whitened'%#ok<NOPTS>
        toc
    end
    clear Filt_hs rez
    
    
    % Load tSp
    tSp = load([dataDir,'\Analysis_April2021\','tSp_',type,'.mat']);
    tSp=tSp.tSp;
    'Cluster and spike data loaded' %#ok<NOPTS>
    toc
    
    window=40;
    AvgSpike = nan(length(peakCh_WMMap),(window*2)+1);
    SEMSpike = nan(length(peakCh_WMMap),(window*2)+1);
    
    for n=1:length(peakCh_WMMap)
        
        NumSpikes=length(tSp{n});
        Spikes = nan(NumSpikes-1,(window*2)+1);
        for s=2:NumSpikes-1
            
            Spikes(s,:)=White_hs(n,(tSp{n}(s)-window):(tSp{n}(s)+window));
            
        end
        AvgSpike(n,:)=mean(Spikes,'omitnan');
        SEMSpike(n,:)=std(Spikes,'omitnan')/sqrt(NumSpikes);
    end
    clear White_hs

    
    
    for m=1:length(peakCh)
        
        subplot(ceil(length(peakCh)/4),4,m)
        
        plot((-window:window),AvgSpike(m,:),'Color',[0 0 0],'LineWidth',2)
        hold on
        fill([-window:(window) fliplr(-window:(window))], [AvgSpike(m,:)-SEMSpike(m,:) fliplr(AvgSpike(m,:)+SEMSpike(m,:))], 'r')
        alpha(0.25)
        title(['Cluster ' int2str(uniqueID{m}) '     Peak Channel ' int2str(peakCh(m)) ], 'Interpreter' , 'none');
        
    end
    waveform = struct;
    waveform.AvgSpike = AvgSpike;
    waveform.SEMSpike = SEMSpike;
    save([saveDir '\Waveform_' type '.mat'],'-struct','waveform')
    saveas(gcf,[saveDir '\waveforms.png'],'png')
    
    close all
    clearvars -except nchans b1 a1 j animalData nSessions type
end