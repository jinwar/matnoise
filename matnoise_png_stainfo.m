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
%array_bgtime = (datenum(2010,12,2) - datenum(1970,1,1))*24*3600;
array_endtime = (datenum(2011,8,1) - datenum(1970,1,1))*24*3600;
time_interval = 3600;
timegrids = array_bgtime:time_interval:array_endtime;
NT = length(timegrids);

for ista = 1:length(stainfo)

	% Initializing the datacover vector
	stainfo(ista).datacover = zeros(length(timegrids)-1,1);

	% Generate sub database that only contain the information of this station
	substr=sprintf('(sta =~/%s/)',stainfo(ista).staname);
	stadbwf=dbsubset(dbwf,substr);
	if dbquery(stadbwf,'dbRECORD_COUNT')<1
		disp(['No data found for station ', stainfo(ista).staname]);
		continue;
	end
	wfbgtime = dbgetv(stadbwf,'time');
	wfendtime = dbgetv(stadbwf,'endtime');

	for itime = 1:NT-1
		ts = timegrids(itime);
		te = timegrids(itime+1);
		% try to find the records that contains this time interval
		segid = find(wfbgtime < ts & wfendtime > te);
		if length(segid) == 1
			stainfo(ista).datacover(itime)=1;
			continue;
		end
		% if not exist, find two segments that connected in this window
		segid = find((wfbgtime > ts & wfbgtime < te) | (wfendtime > ts & wfendtime < te ));
		if length(segid) >1
			bgtimes = sort(wfbgtime(segid));
			endtimes = sort(wfendtime(segid));
			if min(bgtimes) < ts && max(endtimes)>te && (max(bgtimes)-min(endtimes)) < 0.1
				stainfo(ista).datacover(itime)=2;
				continue;
			end
		end
	end % end of time loop
	disp(['Station: ', stainfo(ista).staname , ' Data Coverage:', ...
		num2str(sum(stainfo(ista).datacover)/NT*100),'%' ]);
end % end of station loop

save('stainfo_BHZ','stainfo')
