% This script is used to

clear

% Parameter settings
StableSeg = 1;   % How many hour the instrument stable after being powered
TimeSegN = 1000;

% Initial the dataset
dbpath = 'cdpapuall';
component = 'BHZ';
db = dbopen(dbpath,'r');
dbwf=dblookup_table(db,'wfdisc');
  subsetcomp=sprintf('chan=~/%s/',component);
dbwf=dbsubset(dbwf,subsetcomp);
dbsi=dblookup_table(db,'site');
load stainfo_BHZ.mat

for itime = 1:1
	

	xcor_sum=0;
	xcornum=0;
	for iseg = itime*1000:(itime+1)*1000-1
		if mod(iseg,10)==0
			disp(iseg)
		end
		iscor = 1;
		if stainfo(ista).datacover(iseg) == 0 || stainfo(jsta).datacover(iseg) == 0 
			iscor =0;
		end % Both station have data in this hour
		if sum(stainfo(ista).datacover(iseg-StableSeg+1:iseg)) < StableSeg
			iscor =0;
		end % Make sure station i is stable
		if sum(stainfo(jsta).datacover(iseg-StableSeg+1:iseg)) < StableSeg
			iscor =0;
		end % Make sure station i is stable

		

		% Pre-process the data from two stations
		if iscor
			fftd1 = fft(d1);
			fftd2 = fft(d2);
			% Remove instrument response
			fftd1 = fftd1.*stainfo(ista).resp;
			fftd2 = fftd2.*stainfo(jsta).resp;
			% Whiten
			fftd1 = spectrumwhiten(fftd1);
			fftd2 = spectrumwhiten(fftd2);
			% inverse fft
			d1 = real(ifft(fftd1));
			d2 = real(ifft(fftd2));
			% Down sample to 1s
			d1 = downsample(d1,1/stainfo(ista).dtr);
			d2 = downsample(d2,1/stainfo(ista).dtr);
			if length(d1)~=length(d2)
				disp(['Timeseg: ',num2str(iseg), ...
					'Down Sample Length error' ]);
				iscor=0;
			end
		end 

		% Do Xcor!
		if iscor
			this_xcorr = xcorr(d1,d2);
			this_xcorr=this_xcorr/max(abs(this_xcorr));
			xcor_sum = xcor_sum + this_xcorr;
			xcornum = xcornum+1;
		end
	end % Loop of time segments

	save(sprintf('xcor/%s_%s_%d.mat',stainfo(ista).staname,stainfo(jsta).staname,itime),'xcor_sum','xcornum');
end  % end of large time seg
