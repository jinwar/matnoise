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
for ista = 7:length(stainfo)
	stainfo(ista).staname
	% Number of segments to be averaged.
	avgN = 100;
	N = length(timegrids)-2;
	clear segdata;


	% Generate 100 random segments with data to average
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
	end

	ori_segdata = segdata;
	for iseg = 1:avgN
		segdata(iseg,:) = detrend(segdata(iseg,:),'constant');
		segdata(iseg,:) = detrend(segdata(iseg,:));
		ffty = fft(segdata(iseg,:));
		ffty = ffty.*stainfo(ista).resp';
		segdata(iseg,:) = real(ifft(ffty));
		% Normalize each hour
		if max(abs(segdata(iseg,:))) > 0
			segdata(iseg,:) = segdata(iseg,:)./max(abs(segdata(iseg,:)));
		end
		segdata(iseg,:) = detrend(segdata(iseg,:));
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
	print('-dpng',sprintf('avgspectrum_%s',stainfo(ista).staname));

end  % end of station loop

dbclose(db);
