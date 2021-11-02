% path
% animal
% group
% genotype - wt, LacZ, Cre
% box - 1,2,3,4
% stage
% session
% ChNAcc - usually 0:31
% ChVTA - usually 32:63

function VTA_NAccDataDef

cd('Z:\Ephys Data\postPhy\Anna-Lena\')
group_folders = dir(cd);
f=1;
for g = 3:length(group_folders) % first two entries are \. and \..
    folders = dir([group_folders(g).folder,'\',group_folders(g).name]);
    
    for i = 3:length(folders)
        name = folders(i).name;
        path{f} = [folders(i).folder,'\',name]; %
        animal{f} = name(1:4); %
        group{f} = group_folders(g).name(2); %
        %         genotype{i} = 'wt';
        %         box{i} = 1;
        paradigm = name((length(name)-5):length(name)); %
        if ~isempty(strfind(paradigm,'Test'))
            stage{f} = '3';
            session{f} = name(length(name));
        else
            stage{f} = name(length(name)-3);
            session{f} = name(length(name));
        end       
            ChNAcc{f} = [32:63]; % These values are based on the peak channel labels from Phy and start at 0
            ChVTA{f} = [0:31];

        f=f+1;
    end
end


realInds = find(~cellfun(@isempty,animal));
path=path(realInds);
animal=animal(realInds);
% genotype = genotype(realInds);
% box = box(realInds);
group=group(realInds);
stage=stage(realInds);
session=session(realInds);
ChNAcc=ChNAcc(realInds);
ChVTA=ChVTA(realInds);

save('C:\Users\LabAdmin\Dropbox\ANALYSIS\Anna-Lena\animalDataNAccVTA','path','animal','group','stage','session','ChNAcc','ChVTA'); % 'genotype','box'
end
