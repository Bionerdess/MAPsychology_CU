function trials = getTrials(EVENTS,dataFields,stage,fs)
%%% identify trial events (trial type, CS on, CS off, etc.)

if isfield(EVENTS,'CS_A')
    CS_A = EVENTS.CS_A;
end
if isfield(EVENTS,'CS_B')
    CS_B = EVENTS.CS_B;
end
if isfield(EVENTS,'US_A')
    US_A = EVENTS.US_A;
end
if isfield(EVENTS,'US_B')
    US_B = EVENTS.US_B;
end

A_ITI = [];
ACS_on = [];
AUS_on = [];
ACS_off = [];
B_ITI = [];
BCS_on = [];
BUS_on = [];
BCS_off = [];

for i=1:length(CS_A.ts_start)
    for r=1:length(US_A.ts_start)
       if US_A.ts_start(r)<CS_A.ts_end(i) && US_A.ts_start(r)>CS_A.ts_start(i)
           ITI = CS_A.ts_start(i)-(10*fs);
           A_ITI = [A_ITI;ITI];
           ACS_on = [ACS_on;CS_A.ts_start(i)];
           AUS_on = [AUS_on;US_A.ts_start(r)];
           ACS_off = [ACS_off;CS_A.ts_end(i)];
       end
    end
end
if isfield(EVENTS,'CS_B')
    for i=1:length(CS_B.ts_start)
        for s=1:length(US_B.ts_start)
            if US_B.ts_start(s)<CS_B.ts_end(i) && US_B.ts_start(s)>CS_B.ts_start(i)
                ITI = CS_B.ts_start(i)-(10*fs);
                B_ITI = [B_ITI;ITI];
                BCS_on = [BCS_on;CS_B.ts_start(i)];
                BUS_on = [BUS_on;US_B.ts_start(s)];
                BCS_off = [BCS_off;CS_B.ts_end(i)];
            end
        end
    end
end

for df=1:length(dataFields)
    switch dataFields{df}
        case 'ITI'
            if ~isempty(A_ITI)
                A.(dataFields{df}) = A_ITI;
            end
            if ~isempty(B_ITI)
                B.(dataFields{df}) = B_ITI;
            end
        case 'CS_on'
            if ~isempty(ACS_on)
                A.(dataFields{df}) = ACS_on;
            end
            if ~isempty(BCS_on)
                B.(dataFields{df}) = BCS_on;
            end
        case 'US_on'
            if ~isempty(AUS_on)
                A.(dataFields{df}) = AUS_on;
            end
            if ~isempty(BUS_on)
                B.(dataFields{df}) = BUS_on;
            end
        case 'CS_off'
            if ~isempty(ACS_off)
                A.(dataFields{df}) = ACS_off;
            end
            if ~isempty(BCS_off)
                B.(dataFields{df}) = BCS_off;
            end
    end
end

clear A_ITI ACS_on AUS_on ACS_off B_ITI BCS_on BUS_on BCS_off

if stage=='3'
    for r=1:length(A.CS_on)
        reminderA(r) = find(CS_A.ts_start==A.CS_on(r));
    end
    for s=1:length(B.CS_on)
        reminderB(s) = find(CS_B.ts_start==B.CS_on(s));
    end
    
    A_ITI = [];
    ACS_on = [];
    AUS_omit = [];
    A_CS_off = [];
    B_ITI = [];
    BCS_on = [];
    BUS_omit = [];
    BCS_off = [];
    for i=reminderA(length(reminderA))+1:length(CS_A.ts_start)
        ITI = CS_A.ts_start(i)-(10*fs);
        A_omit = CS_A.ts_start(i)+(9.5*fs);
        A_ITI = [A_ITI;ITI];
        ACS_on = [ACS_on;CS_A.ts_start(i)];
        AUS_omit = [AUS_omit;A_omit];
        A_CS_off = [A_CS_off;CS_A.ts_end(i)];
    end
    for j=reminderB(length(reminderB))+1:length(CS_B.ts_start)
        ITI = CS_B.ts_start(j)-(10*fs);
        B_omit = CS_B.ts_start(j)+(9.5*fs);
        B_ITI = [B_ITI;ITI];
        BCS_on = [BCS_on;CS_B.ts_start(j)];
        BUS_omit = [BUS_omit;B_omit];
        BCS_off = [BCS_off;CS_B.ts_end(j)];
    end
    A_nr.ITI = A_ITI;
    A_nr.CS_on = ACS_on;
    A_nr.US_on = AUS_omit;
    A_nr.CS_off = A_CS_off;
    B_nr.ITI = B_ITI;
    B_nr.CS_on = BCS_on;
    B_nr.US_on = BUS_omit;
    B_nr.CS_off = BCS_off;
end

if exist('A','var')
    trials.A = A;
end
if exist('B','var')
    trials.B = B;
end
if exist('A_nr','var')
    trials.A_nr = A_nr;
end
if exist('B_nr','var')
    trials.B_nr = B_nr;
end
if ~exist('CS_A','var') && ~exist('CS_B','var')
    disp(['ERROR: No rewarded or shocked trials; Check!'])
end