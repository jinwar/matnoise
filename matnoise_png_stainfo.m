%% This code is a modification version of Josh's code "matnoise_v3_batch.m"
%  to do the ambient noise measurement for PNG dataset.
%  Author: Ge Jin, jinwar@gmail.com

% Parameter setting
%
dbpath = 'cdpapuall';
component = 'BHZ';

% Initial the dataset
db = dbopen(dbpath,'r');
dbwf=dblookup_table(db,'wfdisc');
  subsetcomp=sprintf('chan=~/%s/',component);
dbwf=dbsubset(dbwf,subsetcomp);
dbsi=dblookup_table(db,'site');

%gather stn info
stanames = dbgetv(dbsi,'sta');
for ista = 1:length(stanames)
	stan = char(stanames(ista));
	ksite = dbfind(dbsi,sprintf('sta=~/%s/',stan));
	dbsi.record = ksite;
	stainfo(ista).staname = stan;
	stainfo(ista).ksite = ksite;
	stainfo(ista).lat = dbgetv(dbsi,'lat');
	stainfo(ista).lon = dbgetv(dbsi,'lon');
end

% Generate time mark array and mark it for all the stations.
% For PNG dataset, the database time window is:
% 3/02/2010 (061)  6:59:59.30940 - 8/01/2011 (213)  0:00:36.04000
array_bgtime = (datenum(2010,3,2) - datenum(1970,1,1))*24*3600;
array_endtime = (datenum(2011,8,1) - datenum(1970,1,1))*24*3600;
time_interval = 3600;
timegrids = array_bgtime:time_interval:array_endtime;

for ista = 1:length(stainfo)
	stainfo(ista).datacover = zeros(length(timegrids)-1,1);
	for itime = 1:length(timegrids)-1
		% try to find the records that has overlap with this time interval
		ts = timegrids(itime);
		te = timegrids(itime+1);
		substr1=sprintf('(sta =~/%s/) && ((wfdisc.time <= %d && wfdisc.endtime > %d ) || (wfdisc.time > %d && wfdisc.time <%d) || (wfdisc.endtime > %d && wfdisc.endtime <%d) )',stainfo(ista).staname,ts,te,ts,te,ts,te);
		dbtr1=dbsubset(dbwf,substr1);
		if (dbquery(dbtr1,'dbRECORD_COUNT'))==0
			continue;
		end
		% request dataset pointer
		trptr1=trload_css(dbtr1,ts,te);
		% splice segments together
		trsplice(trptr1,50);
		% if the trace number is larger than 1, means the segments cannot be spliced.
		if (dbnrecs(trptr1)~=1)
			continue;
		end
		% Check the start time and end time are correct
		if dbgetv(trptr1,'time')~=ts || dbgetv(trptr1,'endtime')~=te
			continue;
		end
		stainfo(ista).datacover(itime)=1;
	end
end
