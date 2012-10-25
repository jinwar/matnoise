% This script is used to

clear

% Parameter settings
StableHour = 1;   % How many hour the instrument stable after being powered

% Initial the dataset
dbpath = 'cdpapuall';
component = 'BHZ';
db = dbopen(dbpath,'r');
dbwf=dblookup_table(db,'wfdisc');
  subsetcomp=sprintf('chan=~/%s/',component);
dbwf=dbsubset(dbwf,subsetcomp);
dbsi=dblookup_table(db,'site');
load stainfo_BHZ.mat

ista = 1;
jsta = 2;

for iseg = 1:length(timegrids)


end
