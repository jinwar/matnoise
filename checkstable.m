%% Script to check the average spectrum of one stations after being shut down and repowered
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
ista=6;
Nth=3;
stainfo(ista).staname

substr = sprintf('sta=~/%s/',stainfo(ista).staname);
dbtr_site = dbsubset(dbwf,substr);

% Generate 10 random segments with data to average
segnum = 0;
for iseg = Nth:length(stainfo(ista).datacover)
	%if stainfo(ista).datacover(iseg) ==1 && stainfo(ista).datacover(iseg-1) ==0 
	if sum(stainfo(ista).datacover(iseg-Nth+1:iseg)) ==Nth  ...
			&& stainfo(ista).datacover(iseg-Nth) ==0
		if mod(segnum,10)==0
			disp(segnum)
		end
		ts = timegrids(iseg);
		te = timegrids(iseg+1);
		trptr = trload_css(dbtr_site,ts,te);
		trsplice(trptr,5);
		if dbquery(trptr,'dbRECORD_COUNT')~=1
			disp('Cannot merge segments');
			trdestroy(trptr);
			continue;
		end
		d = trextract_data(trptr);
		if length(d)~=180000
			trdestroy(trptr);
			disp('Wrong size of data');
			continue;
		end
		segnum = segnum+1;
		segdata(segnum,:) = d;
		trdestroy(trptr);

	end
end

ori_segdata = segdata;
for iseg = 1:segnum
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
segnum
figure(1)
clf
plot(t,avgsegdata)

figure(2)
clf
loglog(f,abs(fft(avgsegdata)));
xlim([0 10])
xlabel('Frequency (Hz)');

figure(3)
clf
plot(segdata');

dbclose(db);
