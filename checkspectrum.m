%% Script to check the average spectrum of one stations for 100 random hours
%  Written by Ge Jin, jinwar@gmail.com
%  Oct, 2012
clear
% Initial the dataset
dbpath = 'cdpapuall';
component = 'BHZ';
db = dbopen(dbpath,'r');
dbwf=dblookup_table(db,'wfdisc');
  subsetcomp=sprintf('chan=~/%s/',component);
dbwf=dbsubset(dbwf,subsetcomp);
dbsi=dblookup_table(db,'site');
load stainfo_BHZ.mat

% First try to get some normal spectrum
% station AGAN seems to be good (ista=1)
ista=41;
stainfo(ista).staname
% Number of segments to be averaged.
avgN = 100;
N = length(timegrids)-2;


% Generate 10 random segments with data to average
avgind = round(N*rand(avgN,1))+1;
while sum(stainfo(ista).datacover(avgind))<avgN-0.5
	ind = find(stainfo(ista).datacover(avgind)==0);
	avgind(ind) = round(N*rand(length(ind),1))+1;
end
for iseg = 1:avgN
	disp(iseg/avgN)
	substr = sprintf('sta=~/%s/',stainfo(ista).staname);
	dbtr_site = dbsubset(dbwf,substr);
	ts = timegrids(avgind(iseg));
	te = timegrids(avgind(iseg)+1);
	trptr = trload_css(dbtr_site,ts,te);
    trsplice(trptr,5);
	segdata(iseg,:) = trextract_data(trptr);
	trdestroy(trptr);
	% Normalize each hour
	segdata(iseg,:) = detrend(segdata(iseg,:));
	segdata(iseg,:) = segdata(iseg,:)./max(segdata(iseg,:));
end

avgsegdata = mean(segdata,1);
t = 0.02:0.02:3600;
f = [0:1/3600:1/2/0.02,-1/2/0.02+1/3600:1/3600:-1/3600];
figure(1)
clf
plot(t,avgsegdata)

figure(2)
clf
loglog(f,abs(fft(avgsegdata)));
xlim([0 10])
xlabel('Frequency (Hz)');

dbclose(db);
