%% This code is a modification version of Josh's code "matnoise_v3_batch.m"
%  to do the ambient noise measurement for PNG dataset.
%  Author: Ge Jin, jinwar@gmail.com

% Parameter setting
%
dbpath = 'cdpapuall';
component = 'BHZ';

% Initial the dataset

db = dbopen(cdpath,'r');
dbwf=dblookup_table(db,'wfdisc');
  subsetcomp=sprintf('chan=~/%s/',component);
dbwf=dbsubset(dbwf,subsetcomp);
dbsi=dblookup_table(db,'site');

%gather stn info
slats=zeros(1,length(stations)); slons=zeros(1,length(stations));
for ind=1:length(stations)
      thisstn=char(stations(ind));
      ksite=dbfind(dbsi,sprintf('sta=~/%s/',thisstn));
      dbsi.record=ksite;
      [slats(ind) slons(ind)]=dbgetv(dbsi,'lat','lon');
end

