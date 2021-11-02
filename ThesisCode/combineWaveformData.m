%% combine waveform data

clear all

% get data dir
animalData = load('C:\Users\LabAdmin\Dropbox\ANALYSIS\Anna-Lena\animalDataNAccVTA'); %Load the Data Def
nSessions = length(animalData.path);
type = 'NAccVTA_fear';
stage = {'stage1','stage2','stage3'}; % 
saveDir = 'C:\Users\LabAdmin\Dropbox\ANALYSIS\Anna-Lena\NAccVTA_April2021\';

for i=1:length(stage)
    w=struct;
    for s=1:nSessions
        dataDir=[animalData.path{s},'\Analysis_April2021\'];
        cellInfo=load([dataDir,'CellInfo_',type,'.mat']);
        
        switch stage{i}
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
        
        waveform=load([dataDir,'Waveform_',type,'.mat']);
        waveAvg=waveform.AvgSpike;
        waveSEM=waveform.SEMSpike;
        session=['session',cellInfo.session{1}];
        animal = cellInfo.animal{1};
        
        [nCells,~]=size(waveAvg);
        nNeurons=length(cellInfo.CellID);
        
        if nCells ~= nNeurons
            disp('something is wrong! check nCells')
        end
        
        reps=0;
        if isfield(w,session)
            if isfield(w.(session),animal)
                reps=length(w.(session).(animal).AvgSpike);
            end
        end
        
        for c=1:nCells
            w.(session).(animal).AvgSpike{reps+1}{c}=waveAvg(c,:);
            w.(session).(animal).SEMSpike{reps+1}{c}=waveSEM(c,:);
        end
    end
    save([saveDir,'Waveform_',stage{i},'.mat'],'-struct','w')
end