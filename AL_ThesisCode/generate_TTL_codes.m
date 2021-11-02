clear all

TTL_codes = struct;
TTL_codes.TTL_event_codes = [1;2;3;4;5;6;7;8];
% TTL_codes.TTL_codes.TTL_event_names = cell(8,1);
TTL_codes.TTL_event_names{1,1} = 'houseLight';
TTL_codes.TTL_event_names{2,1} = 'US_A';
TTL_codes.TTL_event_names{3,1} = 'CS_A';
TTL_codes.TTL_event_names{4,1} = 'HE';
TTL_codes.TTL_event_names{5,1} = 'CS_B';
TTL_codes.TTL_event_names{6,1} = 'NaN';
TTL_codes.TTL_event_names{7,1} = 'US_B';
TTL_codes.TTL_event_names{8,1} = 'NaN';


cd('C:\Users\LabAdmin\Dropbox\EphysCode\')
save TTL_codes
