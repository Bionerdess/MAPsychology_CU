function EVENTS = getEvents(dataDir)
%%% identify time stamps for events on all TTL channels

cd(dataDir);
% Load TTL map
load('C:\Users\LabAdmin\Dropbox\EphysCode\TTL_codes.mat');
TTL_event_names = TTL_codes.TTL_event_names;

% WM load digital data
file = dir('*Digital*.bin');
fid = fopen(file.name,'r');
t = fread(fid,1,'uint64=>int64'); % timestamp of first sample
din = fread(fid,'int64=>int64'); % digital data

% expand TTLs to binary 
NumChannel = 8; %that correspond to each box - OEP saves each box TTLs to it's 1-8
bit=1;
tic
TTL = false(length(din),NumChannel);
for c=1:NumChannel
    TTL(:,c)= bitget(din,c,'int64');
    bit=bit+1;
end



% Extract list of event names and times

EVENTS = struct;
for r=1:length(TTL_event_names)
    ind = find(TTL(:,r)==0);
    if ~isempty(ind)
        indEvent = find(diff(ind)>10);
        ts_start = [];
        ts_end = [];
        ts_start = ind(1);
        ts_start = [ts_start;ind(indEvent+1)];
        ts_end = ind(indEvent);
        ts_end = [ts_end;ind(length(ind))];
        
        EVENTS.(TTL_event_names{r}).ts_start = ts_start;
        EVENTS.(TTL_event_names{r}).ts_end = ts_end;
    end
    clear ind indEvent
    toc
end

