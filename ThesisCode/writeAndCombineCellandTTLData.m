%% TTL associated neural firing
% Extract animal and cell information as well as timestamps for trial,
% spike and magazine entry events

clear all

% get data dir
animalData = load('C:\Users\LabAdmin\Dropbox\ANALYSIS\Anna-Lena\animalDataNAccVTA'); %Load the Data Def
nSessions = length(animalData.path);
type = 'NAccVTA_fear';

for i = 1:nSessions
    
    dataDir = animalData.path{i};
    cd (dataDir)
    disp(['Processing data of ',dataDir])
    
    saveDir = [dataDir,'\Analysis_April2021'];
    if ~exist(saveDir,'dir')
        mkdir(saveDir)
    end
    tic
    animal = animalData.animal{i};
    animalNum = str2double(animal(length(animal)));
    EVENTS = getEvents(dataDir);
    ts_HE = EVENTS.HE;
    
    'TTL data processed!'
    toc
    
    % Trial Types
    load('C:\Users\LabAdmin\Dropbox\KiloSort\ChannelMap.mat')
    dataFields = {'ITI','CS_on','US_on','CS_off'};
    stage = animalData.stage{i};
    trials = getTrials(EVENTS,dataFields,stage,fs);
    
    toc
    disp('Trials defined')
    
    % Neuron spike timing
    
    [spikeTimes,spikeTimes_bin,CellID,peakCh] = gettSp_bin(dataDir); % get spike timestamps, Cell ID (cluster #), peak Channel
    
    nCells = length(CellID);
    
    disp(['Processed Neuron Data! You got ' int2str(nCells) ' units!'])
    toc
    uniqueID = nan(1,nCells);
    for b=1:nCells
        uniqueID(b) = (CellID(b)*i)+animalNum;
    end
    
    dataToWrite = {'trials','ts_HE','animal','path','session','stage','uniqueID','CellID','peakCh','spikeTimes','spikeTimes_bin'};  % 
    CellInfo = struct;
    
    for d=1:length(dataToWrite)
        switch dataToWrite{d}
            case 'trials'
                fields1 = fields(trials);
                for j=1:length(fields1)
                    fields2 = fields(trials.(fields1{j}));
                    for k=1:length(fields2)
                        T = length(trials.(fields1{j}).(fields2{k}));
                        for t=1:T
                            Trials.(fields1{j}).(fields2{k}){t} = trials.(fields1{j}).(fields2{k})(t);
                        end
                    end
                end
            case 'ts_HE'
                fields1 = fields(ts_HE);
                for j=1:length(fields1)
                    HE = length(ts_HE.(fields1{j}));
                    for h=1:HE
                        tsHE.(fields1{j}){h} = ts_HE.(fields1{j})(h);
                    end
                end
            case 'spikeTimes'
                if ~isempty(spikeTimes)
                    tSp = spikeTimes;
                end
            case 'spikeTimes_bin'
                if ~isempty(spikeTimes_bin)
                    tSp_bin = spikeTimes_bin;
                end
                
            case 'CellID'
                if ~isempty(CellID)
                    for c=1:nCells
                        CellInfo.CellID{c} = CellID(c);
                    end
                end
            case 'uniqueID'
                if ~isempty(CellID)
                    for c=1:nCells
                        CellInfo.uniqueID{c} = uniqueID(c);
                    end
                end
            case 'peakCh'
                if ~isempty(peakCh)
                    for c=1:nCells
                        CellInfo.peakCh{c} = peakCh(c);
                    end
                end
            otherwise
                for c=1:nCells
                    CellInfo.(dataToWrite{d})(c) = animalData.(dataToWrite{d})(i);
                end
        end
    end

    save([saveDir '\Trials_' type '.mat'],'-struct','Trials')
    save([saveDir '\tsHE_' type '.mat'],'-struct','tsHE')
    save([saveDir '\CellInfo_' type '.mat'],'-struct','CellInfo')
    if exist('tSp','var')
        save([saveDir '\tSp_' type],'tSp','-v7.3')
    end
    if exist('tSp_bin','var')
        save([saveDir '\tSpBin_' type],'tSp_bin','-v7.3')
    end
    
    disp(['Data saved for ',dataDir])
    clear Trials tsHE tSp tSp_bin CellInfo
end

%% Combine Data 
% For Stage 3: session1-n=#repetitions of test (session1=test1, session4=test4)

clear all
outputDir = 'C:\Users\LabAdmin\Dropbox\ANALYSIS\Anna-Lena\NAccVTA_April2021\'; %Save location
animalData = load('C:\Users\LabAdmin\Dropbox\ANALYSIS\Anna-Lena\animalDataNAccVTA'); %Load the Data Def
nSessions = length(animalData.path);
dataField = {'tSp','tSpBin','tsHE','Trials','CellInfo'}; % 
paradigms = {'stage1','stage2','stage3'};
combinedData = struct;
type = 'NAccVTA_fear';
for p = 1:length(paradigms)
    for df = 1:length(dataField)
        combinedData = [];
        dataToCombine = 0;
        count=[];
        for i = 1:nSessions % -1
            
            dataDir = animalData.path{i};
            analysisDir = [dataDir,'\Analysis_April2021\'];
            cellInfo = load([analysisDir,'CellInfo_',type]); 
            animal = cellInfo.animal{1};
            field = fields(cellInfo);
            if length(field)>1
                
                switch paradigms{p}
                    case 'stage1'
                        if ~(strcmp(cellInfo.stage,'1'))
                            continue
                        end
                    case 'stage2'
                        if ~(strcmp(cellInfo.stage,'2'))
                            continue
                        end
                    case 'stage3'
                        if ~(strcmp(cellInfo.stage,'3'))
                            continue
                        end
                end
                
                sessionName = (['session',animalData.session{i}]);
                [combinedData] = combineMultiSessionData(dataField{df},analysisDir,type,combinedData,sessionName,animal);
                dataToCombine = 1;
            end
            clear field
        end
        
        if dataToCombine
            if ~isstruct(combinedData)% || ismember(dataField{df},{'spikeTimes'})
                data = combinedData;
                saveData = struct;
                saveData.(dataField{df}) = data;
                save([outputDir,dataField{df},'_',paradigms{p}],'-struct','saveData');
            else
                save([outputDir,dataField{df},'_',paradigms{p}],'-struct','combinedData','-v7.3');
            end
        end
        disp([paradigms{p},': ',dataField{df},' Data Combined'])
        
        fclose('all');
        disp('Data Saved')
        
    end
end