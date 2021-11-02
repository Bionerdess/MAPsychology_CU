%% calculate repeated measures ANOVA - select .csv file organized in columns = within factor levels by nfactors (e.g. 2 factors with 8 levels each needs to be 16 columns)
% debug for more than 2 factors!

start_path = 'C:\Users\LabAdmin\Dropbox\ANALYSIS\';
dialog_title = 'Select .csv file containing data';
[filename,pathname] = uigetfile('*.csv','Select data file');

data =importdata([pathname,filename]);
data = data.data;
f = inputdlg({'Enter number of factors','Enter levels of factors'},'Define factors');
f = str2double(f);

[MainA,MainB,Interaction,SimpleAatB,SimpleBatA] = mystats(data,f(1),[f(1) f(2)]);