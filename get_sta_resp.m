%% This function is used to get the station response with the required data length
% Written by Ge Jin
%

dbpath = 'cdpapuall';
time_interval = 6000;
lo_corner =  0.005;
component = 'BHZ';

db = dbopen(dbpath,'r');
dbwf=dblookup_table(db,'wfdisc');
  subsetcomp=sprintf('chan=~/%s/',component);
dbwf=dbsubset(dbwf,subsetcomp);
dbsi=dblookup_table(db,'site');

%gather stn info
stanames = dbgetv(dbsi,'sta');
for ista = 1:length(stanames)
	stan = char(stanames(ista));
%	if ista == 16
%		disp('Change station MAPA to MAPM');
%		stan = 'MAPM';
%	end
	ksite = dbfind(dbsi,sprintf('sta=~/%s/',stan));
	dbsi.record = ksite;
	staresp(ista).staname = stan;
	staresp(ista).ksite = ksite;
	staresp(ista).lat = dbgetv(dbsi,'lat');
	staresp(ista).lon = dbgetv(dbsi,'lon');
end

%gather response
dbsen=dblookup_table(db,'sensor');
dbin=dblookup_table(db,'instrument');
dbsnin=dbjoin(dbsen,dbin);
for ind=1:length(staresp)
	disp(['Getting Instrument Resp for station:' staresp(ind).staname])
	thisstn=char(staresp(ind).staname); 
	clear resp dtr
	[resp dtr]=calcinstresp(dbsnin,thisstn,component,-1, time_interval, lo_corner);
	staresp(ind).resp = resp;
	staresp(ind).dtr = dtr;
end

save staresp.mat staresp
