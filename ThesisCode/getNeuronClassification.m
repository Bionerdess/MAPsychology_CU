%% Neuron identification (dopamine vs. non-dopamine)

clear all

% get data dir
dataDir = 'C:\Users\LabAdmin\Dropbox\ANALYSIS\Anna-Lena\NAccVTA_April2021\'; % combined Data location

groups = {'stage3'};
region = {'VTA','NAcc'};

for g=1:length(groups)
    fh1 = figure;
    fh2 = figure;
    
    ampRatio_VTA = [];
    ampRatio_NAcc = [];
    halfDur_VTA = [];
    halfDur_NAcc = [];
    uniqueVTA = [];
    uniqueNAcc = [];

    waveData = load([dataDir,'Waveform_',groups{g},'.mat']); %load trial data
    CellInfo = load([dataDir,'CellInfo_',groups{g},'.mat']); % load CellInfo data (animal,CellID,path,peakCh,session,stage)
    
    session = fields(CellInfo);
    for s=1:length(session)
        animal = fields(CellInfo.(session{s}));
        
        for i = 1:length(animal)
            reps=length(CellInfo.(session{s}).(animal{i}).uniqueID);
            for j=1:reps
                peakCh = CellInfo.(session{s}).(animal{i}).peakCh{j};
                ID = CellInfo.(session{s}).(animal{i}).uniqueID{j};
                nCells=length(ID);
                
                ampRatio = nan(2,nCells);
                halfDur = nan(2,nCells);
                ID_VTA = nan(1,nCells);
                ID_NAcc = nan(1,nCells);
                
                for c=1:nCells
                    wave = waveData.(session{s}).(animal{i}).AvgSpike{j}{c};
                    points = length(wave);
                    
                    window = (points-61)/2;
                    wave= wave(:,window:points-window);
                    
                    minimum = min(wave); % get absolute min of waveform
                    indMin = find(wave == minimum); % get column (x value) of minimum
                    indMin = indMin(:,ceil(end/2));
                    max1 = max(wave(1:indMin)); % get first max of waveform
                    indMax2 = find(wave == max(wave(indMin:length(wave)))); % get column (x value) of second max
                    indMax2 = indMax2(:,ceil(end/2));
                    if peakCh{c} <32 % if VTA neuron
                        ampRatio(1,c) = (abs(minimum)-abs(max1))/(abs(minimum)+abs(max1));
                        halfDur(1,c) = (indMax2-indMin)/25;
                        ID_VTA(c) = ID{c};
                    elseif peakCh{c} >31 % if NAcc neuron
                        ampRatio(2,c) = (abs(minimum)-abs(max1))/(abs(minimum)+abs(max1));
                        halfDur(2,c) = (indMax2-indMin)/25;
                        ID_NAcc(c) = ID{c};
                        
                    end
                    
                end
                
                ampRatio_VTA = [ampRatio_VTA,rmmissing(ampRatio(1,:))];
                ampRatio_NAcc = [ampRatio_NAcc,rmmissing(ampRatio(2,:))];
                halfDur_VTA = [halfDur_VTA,rmmissing(halfDur(1,:))];
                halfDur_NAcc = [halfDur_NAcc,rmmissing(halfDur(2,:))];
                
                uniqueVTA = [uniqueVTA,rmmissing(ID_VTA)];
                uniqueNAcc = [uniqueNAcc,rmmissing(ID_NAcc)];
            end
        end
    end
    
    VTA_values = cat(1,ampRatio_VTA,halfDur_VTA);
    VTA_values = VTA_values';
    
    [idx_VTA,C_VTA,sumD_VTA,D_VTA] = kmeans(VTA_values,2,'Replicates',200);
    
    VTA_values(:,3) = uniqueVTA'; % third column of VTA values is the cell's unique ID
    V = unique(VTA_values(:,3),'stable');
    
    NAcc_values = cat(1,ampRatio_NAcc,halfDur_NAcc);
    NAcc_values = NAcc_values';
    [idx_NAcc,C_NAcc,sumD_NAcc,D_NAcc] = kmeans(NAcc_values,2,'Replicates',200);
    
    NAcc_values(:,3) = uniqueNAcc'; % third column of NAcc values is the cell's unique ID
    N = unique(NAcc_values(:,3),'stable');
    
    VTA_values = num2cell(VTA_values);
    NAcc_values = num2cell(NAcc_values);
    % get indices of rows belonging to each cluster
    idxV1 = find(idx_VTA ==1);
    idxV2 = find(idx_VTA ==2);
    idxN1 = find(idx_NAcc ==1);
    idxN2 = find(idx_NAcc ==2);
    
    % get distance values for each cluster
    DV1 = D_VTA(idxV1,1);
    DV2 = D_VTA(idxV2,2);
    DN1 = D_NAcc(idxN1,1);
    DN2 = D_NAcc(idxN2,2);
    
    % calc SD for each cluster (D1, D2)
    DV1_sd = std(DV1);
    DV2_sd = std(DV2);
    DN1_sd = std(DN1);
    DN2_sd = std(DN2);
    
    % third column of distance vector will contain SD values
    D_VTA(idxV1,3) = DV1_sd;
    D_VTA(idxV2,3) = DV2_sd;
    D_NAcc(idxN1,3) = DN1_sd;
    D_NAcc(idxN2,3) = DN2_sd;
    
    % standardize all values/distances of all clusters
    D_VTA(:,4) = D_VTA(:,1)./D_VTA(:,3); % fourth column will contain all distances standardized to cluster 1
    D_VTA(:,5) = D_VTA(:,2)./D_VTA(:,3); % fifth column will contain all distances standardized to cluster 2
    D_NAcc(:,4) = D_NAcc(:,1)./D_NAcc(:,3);
    D_NAcc(:,5) = D_NAcc(:,2)./D_NAcc(:,3);
    
    for dv = 1:length(D_VTA)
        if D_VTA(dv,4) <2 && D_VTA(dv,5) <2 % less than 2 SD away from both
            D_VTA(dv,6) = idx_VTA(dv);
            VTA_values{dv,4} = idx_VTA(dv);
        elseif D_VTA(dv,4) <2 && D_VTA(dv,5) >2 % less than 2SD away from cl1 but greater 2SD from cl2 -> belongs to cluster 1
            D_VTA(dv,6) = 1;
            VTA_values{dv,4} = 1;
        elseif D_VTA(dv,5) <2 && D_VTA(dv,4) >2 % less than 2SD away from cl2 but greater 2SD from cl1 -> belongs to cluster 2
            D_VTA(dv,6) = 2;
            VTA_values{dv,4} = 2;
        else
            D_VTA(dv,6) = 3; % greater 2SD away from both
            VTA_values{dv,4} = 3;
        end
    end
    
    for dn = 1:length(D_NAcc)
        if D_NAcc(dn,4) <2 && D_NAcc(dn,5) <2 % less than 2 SD away from both
            D_NAcc(dn,6) = idx_NAcc(dn);
            NAcc_values{dn,4} = idx_NAcc(dn);
        elseif D_NAcc(dn,4) <2 && D_NAcc(dn,5) >2 % less than 2SD away from cl1 but greater 2SD from cl2 -> belongs to cluster 1
            D_NAcc(dn,6) = 1;
            NAcc_values{dn,4} = 1;
        elseif D_NAcc(dn,5) <2 && D_NAcc(dn,4) >2 % less than 2SD away from cl2 but greater 2SD from cl1 -> belongs to cluster 2
            D_NAcc(dn,6) = 2;
            NAcc_values{dn,4} = 2;
        else
            D_NAcc(dn,6) = 3; % greater 2SD away from both
            NAcc_values{dn,4} = 3;
        end
    end
    
    ind1 = find([VTA_values{:,4}]==1);
    ind2 = find([VTA_values{:,4}]==2);
    [h1,p1,ci1,stats1] = ttest2([VTA_values{ind1,1}],[VTA_values{ind2,1}]); % ampRatio
    [h2,p2,ci2,stats2] = ttest2([VTA_values{ind1,2}],[VTA_values{ind2,2}]); % halfDur
    
    % dopaminergic neurons: smaller ampRatio & longer halfDur
    if p1<0.05 && stats1.tstat <0 && p2<0.05 && stats2.tstat>0 % ampRatio cluster1 sig smaller than cluster2 & halfDur cluster1 sig larger than cluster2
        for v=1:length(ind1)
            VTA_values{ind1(v),5} = {'d'};
        end
        for v=1:length(ind2)
            VTA_values{ind2(v),5} = {'nd'};
        end
    elseif p1<0.05 && stats1.tstat>0 && p2<0.05 && stats2.tstat<0 % ampRatio cluster1 sig larger than cluster2 & halfDur cluster1 sig smaller than cluster2
        for v=1:length(ind1)
            VTA_values{ind1(v),5} = {'nd'};
        end
        for v=1:length(ind2)
            VTA_values{ind2(v),5} = {'d'};
        end
    end
    
    figure(fh1)
    plot([VTA_values{D_VTA(:,6)==1,1}],[VTA_values{D_VTA(:,6)==1,2}],'r.','MarkerSize',30)
    hold on
    plot([VTA_values{D_VTA(:,6)==2,1}],[VTA_values{D_VTA(:,6)==2,2}],'b.','MarkerSize',30)
    hold on
    if ~isempty(find(D_VTA(:,6)==0,1))
        plot([VTA_values{D_VTA(:,6)==0,1}],[VTA_values{D_VTA(:,6)==0,2}],'go')
        hold on
    end
    if ~isempty(find(D_VTA(:,6)==3,1))
        plot([VTA_values{D_VTA(:,6)==3,1}],[VTA_values{D_VTA(:,6)==3,2}],'ko')
        hold on
    end
    plot(C_VTA(:,1),C_VTA(:,2),'kx')
    hold on
    
    figure(fh2)
    plot([NAcc_values{D_NAcc(:,6)==1,1}],[NAcc_values{D_NAcc(:,6)==1,2}],'r.','MarkerSize',30)
    hold on
    plot([NAcc_values{D_NAcc(:,6)==2,1}],[NAcc_values{D_NAcc(:,6)==2,2}],'b.','MarkerSize',30)
    hold on
    if ~isempty(find(D_NAcc(:,6)==0,1))
        plot([NAcc_values{D_NAcc(:,6)==0,1}],[NAcc_values{D_NAcc(:,6)==0,2}],'go')
        hold on
    end
    if ~isempty(find(D_NAcc(:,6)==3,1))
        plot([NAcc_values{D_NAcc(:,6)==3,1}],[NAcc_values{D_NAcc(:,6)==3,2}],'ko')
        hold on
    end
    plot(C_NAcc(:,1),C_NAcc(:,2),'kx')
    hold on
    
    figure(fh1)
    title([groups{g} ' - VTA - Dopamine vs. Non-Dopamine'])
    xlabel('Amplitude Ratio [(fn-fp)/(fn+fp)]')
    ylabel('Half Duration [ms]')
    legend({'Cluster 1','Cluster 2','Overlap (<2SD from more than one cluster)','>2SD from a cluster','Cluster centroid'},'Location','northeast')
    saveas(fh1,[dataDir,'images\VTA_Classification_D-nD_',groups{g},'.svg'],'svg')
    
    figure(fh2)
    title([groups{g} ' - NAcc - Dopamine vs. Non-Dopamine'])
    xlabel('Amplitude Ratio [(fn-fp)/(fn+fp)]')
    ylabel('Half Duration [ms]')
    legend({'Cluster 1','Cluster 2','Overlap (<2SD from more than one cluster)','>2SD from a cluster','Cluster centroid'},'Location','northeast')
    saveas(fh2,[dataDir,'images\NAcc_Classification_D-nD_',groups{g},'.svg'],'svg')
    
    save([dataDir,'VTA_values_noOverlap_',groups{g}],'VTA_values')
    save([dataDir,'NAcc_values_noOVerlap_',groups{g}],'NAcc_values')
    
    clear VTA_values NAcc_values
end

