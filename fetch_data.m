function fetch_data(ista,itime)
% Function to fetch 

% Parameter settings
StableSeg = 1;   % How many hour the instrument stable after being powered
TimeSegN = 1000; % Define how large is one mat file segment

% Initial the dataset
dbpath = 'cdpapuall';
component = 'BDH';
db = dbopen(dbpath,'r');
dbwf=dblookup_table(db,'wfdisc');
  subsetcomp=sprintf('chan=~/%s/',component);
dbwf=dbsubset(dbwf,subsetcomp);
dbsi=dblookup_table(db,'site');
load stainfo_BHZ.mat

segdata = zeros(TimeSegN,3600);
for iseg = itime*TimeSegN:(itime+1)*TimeSegN-1;
	if mod(iseg,100)==0
		disp(iseg)
	end
	if iseg == 0 
		continue;
	end
	if stainfo(ista).datacover(iseg) ==1
		substr = sprintf('sta=~/%s/',stainfo(ista).staname);
		dbtr_site = dbsubset(dbwf,substr);
		ts = timegrids(iseg);
		te = timegrids(iseg+1);
		trptr = trload_css(dbtr_site,ts,te);
		trsplice(trptr,5);
		if dbquery(trptr,'dbRECORD_COUNT')~=1
			disp(['Timeseg: ',num2str(iseg),', Station: ',stainfo(ista).staname,...
				'Cannot merge segments' ]);
			trdestroy(trptr);
			continue;
		end
		d1 = trextract_data(trptr);
		trdestroy(trptr);
		if length(d1) ~=180000
			disp(['Timeseg: ',num2str(iseg),', Station: ',stainfo(ista).staname,...
				'Data Length error' ]);
			continue;
		end

		d1 = downsample(d1,1/stainfo(ista).dtr);
		if length(d1)~=3600
			disp(['Timeseg: ',num2str(iseg),', Station: ',stainfo(ista).staname,...
				'Downsample Data Length error' ]);
			continue;
		end
		segdata(iseg,:)=d1;
	end
end % Loop of time segments

save(sprintf('data/%s_%d_%s.mat',stainfo(ista).staname,itime,component),'segdata');
