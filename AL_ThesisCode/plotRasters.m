%% plot raster plots

clear all

% get data dir
dataDir = 'C:\Users\LabAdmin\Dropbox\ANALYSIS\Anna-Lena\NAccVTA_April2021\'; % combined Data location

groups = {'stage3'}; 
region = {'VTA','NAcc'}; %
type = {'A_nr','B_nr'};

for g=1:length(groups)
    saveDir = [dataDir,'images\RasterPlots'];
    trial_data = load([dataDir,'Trials_',groups{g}]); %load trial data
    CellInfo = load([dataDir,'CellInfo_',groups{g}]); % load CellInfo data (animal,CellID,path,peakCh,session,stage)
    tSp = load([dataDir,'tSp_',groups{g}]); % load spike timestamps
    session = fields(trial_data);
    for s=1:length(session)
        for f=1:length(type)
            animal = fields(trial_data.(session{s}).(type{f}));
            for i = 1:length(animal)
                trial = fields(trial_data.(session{s}).(type{f}).(animal{i}));
                reps=length(trial_data.(session{s}).(type{f}).(animal{i}).(trial{1}));
                for j=1:reps
                    events = trial_data.(session{s}).(type{f}).(animal{i});
                    nTrials = length([events.(trial{1}){j}]);
                    event_ts = nan(nTrials,length(trial));
                    for e=1:length(trial)
                        event_ts(:,e) = cell2mat(events.(trial{e}){j})-cell2mat(events.(trial{1}){j});
                    end
                    nCells = length([tSp.(session{s}).(animal{i}){j}]);
                    peakCh = CellInfo.(session{s}).(animal{i}).peakCh{j};
                    uniqueID = CellInfo.(session{s}).(animal{i}).uniqueID{j};
                    for r=1:length(region)
%                         count=1;
                        for c=1:nCells
                            
                            switch region{r}
                                case 'VTA'
                                    if peakCh{c}>31 % If NAcc peakCh (in phy file from 32-63) continue to next switch
                                        continue
                                    end
                                case 'NAcc'
                                    if peakCh{c}<32 % If VTA peakCh (in phy file from 0-31) continue to next switch
                                        continue
                                    end
                            end
                            
                            fh1 = figure;
                            spikes = [tSp.(session{s}).(animal{i}){j}{c}];
                            PSTH = zeros(nTrials,200); % 20sec trial*100ms windows = 2000 bins
                            for t=1:nTrials
                                start = events.ITI{j}{t}+5*25000;
                                stop = events.CS_off{j}{t}+5*25000;
                                ind = find(spikes>start & spikes<stop);
                                data = spikes(ind)-start;
                                for x=0:199 % for 100ms bins (100ms = 2500 timebins); for 1 sec bins (1 sec = 25000 timestamps)
                                    PSTH(t,x+1) = length(find(data>= x*2500 & data< x*2500+2500));
                                end
                                
                                if length(data)>2
                                    figure(fh1)
                                    subplot(2,1,2)
                                    line([data data],[t-0.4 t+0.4],'Color','k','LineWidth',0.02); %
                                    hold on
                                end
                            end
                            PSTH=sum(PSTH,1);
                            y=[50 145 150];
                            figure(fh1)
                            subplot(2,1,1)
                            bar(PSTH,1)
                            hold on
                            for w=1:length(y)
                            line([y(w) y(w)],[0 max(PSTH)+1],'Color','r','LineStyle','--','LineWidth',1)
                            end
                            
                            subplot(2,1,2)
                            for e=2:length(event_ts(1,:))
                                line([event_ts(1,e)-5*25000 event_ts(1,e)-5*25000],[0 nTrials+0.5],'Color','r','LineWidth',1.5)
                                hold on
                            end
                            
                            figure(fh1)
                            subplot(2,1,1)
                            axis([0 200 0 max(PSTH)+2])
                            ylabel('nSpikes')
                            title(['Cell ',num2str(uniqueID{c})])
                            
                            subplot(2,1,2)
                            set(gca, 'YDir','reverse')
                            set(gca,'xtick',[])
                            set(gca,'box','on')
                            axes(gca)
                            title([region{r},' (' num2str(peakCh{c}) ') ' type{f} ': ' groups{g} ' ' session{s}])
                            xlabel('    ITI       CS on                  US     CS off')
                            ylabel('Trial')
                            
                            saveas(fh1,[saveDir,'\',region{r},animal{i},'se',num2str(s),'Cell',num2str(uniqueID{c}),'_',type{f},'.svg'],'svg');
                            close all
                        end
                    end
                    
                end
            end
        end
    end
end
