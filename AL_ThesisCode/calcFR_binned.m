%% calc Firing Rate in 100 ms bins (by trial)

clear all

% get data dir
dataDir = 'C:\Users\LabAdmin\Dropbox\ANALYSIS\Anna-Lena\NAccVTA_April2021\'; % combined Data location

groups = {'stage3'};
% region = {'VTA','NAcc'};
bin_size = 100; % bin size in ms
trial_length = 35000; % in ms = seconds*1000 (ITI + baseline + CS + US/post-CS)
bin_width = 25000*(bin_size/1000); % number of samples during one bin (sampling rate = 25000Hz * 0.05 sec bin_size)
% 
for g=1:length(groups)
    FR = struct;
    VTA_IDs = [];
    NAcc_IDs = [];
    trials = load([dataDir,'Trials_',groups{g}]); %load trial data
    CellInfo = load([dataDir,'CellInfo_',groups{g}]); % load CellInfo data (animal,CellID,path,peakCh,session,stage)
    tSp = load([dataDir,'tSp_',groups{g}]); % load spike timestamps
    session = fields(trials);
    count = 1;
    for s=1:length(session)
        type = fields(trials.(session{s})); % organize each session by Trials x events (trial_start,CS_onset,US,CS_offset)
        for f=1:length(type)
            animal = fields(trials.(session{s}).(type{f}));

            for i = 1:length(animal)
                trial = fields(trials.(session{s}).(type{f}).(animal{i}));
                reps=length(trials.(session{s}).(type{f}).(animal{i}).(trial{1}));
                for j=1:reps
                    events = trials.(session{s}).(type{f}).(animal{i});
                    nTrials = length([events.(trial{1}){j}]);
                    nCells = length([tSp.(session{s}).(animal{i}){j}]);
                    peakCh = CellInfo.(session{s}).(animal{i}).peakCh{j};
                    uniqueID = CellInfo.(session{s}).(animal{i}).uniqueID{j};
                    
                    for c=1:nCells
                        spikes = [tSp.(session{s}).(animal{i}){j}{c}];
                        data = zeros(nTrials,trial_length/bin_size);
                        
                        for t=1:nTrials
                            start = events.ITI{j}{t}-10*25000;
                            stop = events.CS_off{j}{t}+5*25000;
                            data_ind = find(spikes>=start & spikes<=stop);
                            data_raw = spikes(data_ind)-start;
                            
                            for x=0:size(data,2)-1
                                data(t,x+1) = length(find(data_raw>= x*bin_width & data_raw < x*bin_width+bin_width));
                            end
                            
                        end
                        
                        ITI = data(:,1:10000/bin_size);
                        baseline = data(:,10000/bin_size+1:20000/bin_size);
                        CS = data(:,20000/bin_size+1:30000/bin_size);
                        US = data(:,30000/bin_size+1:trial_length/bin_size);
                        ITI_binavg = mean(ITI,1);
                        ITI_mean = mean(ITI_binavg);
                        ITI_sd = std(ITI_binavg);
                        
                        baseline_z = (baseline-ITI_mean)/ITI_sd;
                        avgBaseline = mean(baseline_z,1);
                        avgBaseline = mean(avgBaseline);
                        CS_z = (CS-ITI_mean)/ITI_sd;
                        US_z = (US-ITI_mean)/ITI_sd;
                        diffScoreCS = CS_z-avgBaseline;
                        diffScoreUS = US_z-avgBaseline;
                        CSraw = CS/(bin_size/1000);
                        USraw = US/(bin_size/1000);
                        BLraw = baseline/(bin_size/1000);
                        
                        All_values{count,1} = uniqueID{c};
                        All_values{count,3} = baseline_z;
                        All_values{count,4} = CS_z;
                        All_values{count,5} = US_z;
                        All_values{count,6} = diffScoreCS;
                        All_values{count,7} = diffScoreUS;
                        All_values{count,8} = CSraw;
                        All_values{count,9} = USraw;
                        All_values{count,10} = BLraw;
                        All_values{count,11} = {type{f}};
                        
                        if peakCh{c}<32 % If VTA cell
                            All_values{count,2} = {'VTA'};
                        end
                        if peakCh{c}>31 % If NAcc cell
                            All_values{count,2} = {'NAcc'};
                        end
                           
                        count = count+1;
                    end
                    
                end
            end
        end
    end
    
    VTA_A_ind = find(cellfun(@(All_values) strcmp(All_values,'VTA'),All_values(:,2))==1 & cellfun(@(All_values) strcmp(All_values,'A'),All_values(:,11))==1);
    VTA_B_ind = find(cellfun(@(All_values) strcmp(All_values,'VTA'),All_values(:,2))==1 & cellfun(@(All_values) strcmp(All_values,'B'),All_values(:,11))==1);
    NAcc_A_ind = find(cellfun(@(All_values) strcmp(All_values,'NAcc'),All_values(:,2))==1 & cellfun(@(All_values) strcmp(All_values,'A'),All_values(:,11))==1);
    NAcc_B_ind = find(cellfun(@(All_values) strcmp(All_values,'NAcc'),All_values(:,2))==1 & cellfun(@(All_values) strcmp(All_values,'B'),All_values(:,11))==1);
    VTA_A_nr_ind = find(cellfun(@(All_values) strcmp(All_values,'VTA'),All_values(:,2))==1 & cellfun(@(All_values) strcmp(All_values,'A_nr'),All_values(:,11))==1);
    NAcc_A_nr_ind = find(cellfun(@(All_values) strcmp(All_values,'NAcc'),All_values(:,2))==1 & cellfun(@(All_values) strcmp(All_values,'A_nr'),All_values(:,11))==1);
    VTA_B_nr_ind = find(cellfun(@(All_values) strcmp(All_values,'VTA'),All_values(:,2))==1 & cellfun(@(All_values) strcmp(All_values,'B_nr'),All_values(:,11))==1);
    NAcc_B_nr_ind = find(cellfun(@(All_values) strcmp(All_values,'NAcc'),All_values(:,2))==1 & cellfun(@(All_values) strcmp(All_values,'B_nr'),All_values(:,11))==1);
    
    
    if ~isempty(VTA_A_ind)
        FR.A.VTA.uniqueID = All_values(VTA_A_ind,1);
        FR.A.VTA.baseline_z = All_values(VTA_A_ind,3);
        FR.A.VTA.CS_z = All_values(VTA_A_ind,4);
        FR.A.VTA.US_z = All_values(VTA_A_ind,5);
        FR.A.VTA.diffScoreCS = All_values(VTA_A_ind,6);
        FR.A.VTA.diffScoreUS = All_values(VTA_A_ind,7);
        FR.A.VTA.CS = All_values(VTA_A_ind,8);
        FR.A.VTA.US = All_values(VTA_A_ind,9);
        FR.A.VTA.baseline = All_values(VTA_A_ind,10);
    end
    if ~isempty(VTA_B_ind)
        FR.B.VTA.uniqueID = All_values(VTA_B_ind,1);
        FR.B.VTA.baseline_z = All_values(VTA_B_ind,3);
        FR.B.VTA.CS_z = All_values(VTA_B_ind,4);
        FR.B.VTA.US_z = All_values(VTA_B_ind,5);
        FR.B.VTA.diffScoreCS = All_values(VTA_B_ind,6);
        FR.B.VTA.diffScoreUS = All_values(VTA_B_ind,7);
        FR.B.VTA.CS = All_values(VTA_B_ind,8);
        FR.B.VTA.US = All_values(VTA_B_ind,9);
        FR.B.VTA.baseline = All_values(VTA_B_ind,10);
    end
    if ~isempty(NAcc_A_ind)
        FR.A.NAcc.uniqueID = All_values(NAcc_A_ind,1);
        FR.A.NAcc.baseline_z = All_values(NAcc_A_ind,3);
        FR.A.NAcc.CS_z = All_values(NAcc_A_ind,4);
        FR.A.NAcc.US_z = All_values(NAcc_A_ind,5);
        FR.A.NAcc.diffScoreCS = All_values(NAcc_A_ind,6);
        FR.A.NAcc.diffScoreUS = All_values(NAcc_A_ind,7);
        FR.A.NAcc.CS = All_values(NAcc_A_ind,8);
        FR.A.NAcc.US = All_values(NAcc_A_ind,9);
        FR.A.NAcc.baseline = All_values(NAcc_A_ind,10);
    end
    if ~isempty(NAcc_B_ind)
        FR.B.NAcc.uniqueID = All_values(NAcc_B_ind,1);
        FR.B.NAcc.baseline_z = All_values(NAcc_B_ind,3);
        FR.B.NAcc.CS_z = All_values(NAcc_B_ind,4);
        FR.B.NAcc.US_z = All_values(NAcc_B_ind,5);
        FR.B.NAcc.diffScoreCS = All_values(NAcc_B_ind,6);
        FR.B.NAcc.diffScoreUS = All_values(NAcc_B_ind,7);
        FR.B.NAcc.CS = All_values(NAcc_B_ind,8);
        FR.B.NAcc.US = All_values(NAcc_B_ind,9);
        FR.B.NAcc.baseline = All_values(NAcc_B_ind,10);
    end
    if ~isempty(VTA_A_nr_ind)
        FR.A_nr.VTA.uniqueID = All_values(VTA_A_nr_ind,1);
        FR.A_nr.VTA.baseline_z = All_values(VTA_A_nr_ind,3);
        FR.A_nr.VTA.CS_z = All_values(VTA_A_nr_ind,4);
        FR.A_nr.VTA.US_z = All_values(VTA_A_nr_ind,5);
        FR.A_nr.VTA.diffScoreCS = All_values(VTA_A_nr_ind,6);
        FR.A_nr.VTA.diffScoreUS = All_values(VTA_A_nr_ind,7);
        FR.A_nr.VTA.CS = All_values(VTA_A_nr_ind,8);
        FR.A_nr.VTA.US = All_values(VTA_A_nr_ind,9);
        FR.A_nr.VTA.baseline = All_values(VTA_A_nr_ind,10);
    end
    if ~isempty(NAcc_A_nr_ind)
        FR.A_nr.NAcc.uniqueID = All_values(NAcc_A_nr_ind,1);
        FR.A_nr.NAcc.baseline_z = All_values(NAcc_A_nr_ind,3);
        FR.A_nr.NAcc.CS_z = All_values(NAcc_A_nr_ind,4);
        FR.A_nr.NAcc.US_z = All_values(NAcc_A_nr_ind,5);
        FR.A_nr.NAcc.diffScoreCS = All_values(NAcc_A_nr_ind,6);
        FR.A_nr.NAcc.diffScoreUS = All_values(NAcc_A_nr_ind,7);
        FR.A_nr.NAcc.CS = All_values(NAcc_A_nr_ind,8);
        FR.A_nr.NAcc.US = All_values(NAcc_A_nr_ind,9);
        FR.A_nr.NAcc.baseline = All_values(NAcc_A_nr_ind,10);
    end
    if ~isempty(VTA_B_nr_ind)
        FR.B_nr.VTA.uniqueID = All_values(VTA_B_nr_ind,1);
        FR.B_nr.VTA.baseline_z = All_values(VTA_B_nr_ind,3);
        FR.B_nr.VTA.CS_z = All_values(VTA_B_nr_ind,4);
        FR.B_nr.VTA.US_z = All_values(VTA_B_nr_ind,5);
        FR.B_nr.VTA.diffScoreCS = All_values(VTA_B_nr_ind,6);
        FR.B_nr.VTA.diffScoreUS = All_values(VTA_B_nr_ind,7);
        FR.B_nr.VTA.CS = All_values(VTA_B_nr_ind,8);
        FR.B_nr.VTA.US = All_values(VTA_B_nr_ind,9);
        FR.B_nr.VTA.baseline = All_values(VTA_B_nr_ind,10);
    end
    if ~isempty(NAcc_B_nr_ind)
        FR.B_nr.NAcc.uniqueID = All_values(NAcc_B_nr_ind,1);
        FR.B_nr.NAcc.baseline_z = All_values(NAcc_B_nr_ind,3);
        FR.B_nr.NAcc.CS_z = All_values(NAcc_B_nr_ind,4);
        FR.B_nr.NAcc.US_z = All_values(NAcc_B_nr_ind,5);
        FR.B_nr.NAcc.diffScoreCS = All_values(NAcc_B_nr_ind,6);
        FR.B_nr.NAcc.diffScoreUS = All_values(NAcc_B_nr_ind,7);
        FR.B_nr.NAcc.CS = All_values(NAcc_B_nr_ind,8);
        FR.B_nr.NAcc.US = All_values(NAcc_B_nr_ind,9);
        FR.B_nr.NAcc.baseline = All_values(NAcc_B_nr_ind,10);
    end
    
    save([dataDir,'FR_',groups{g}],'-struct','FR')
end