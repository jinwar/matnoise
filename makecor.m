% This script is used to

clear

% Parameter settings
StableSeg = 1;   % How many hour the instrument stable after being powered
TimeSegN = 1000;

load stainfo_BHZ.mat

for ista = 1:length(stainfo)
	for jsta = ista+1:length(stainfo)
		for itime = 1:11
			disp([stainfo(ista).staname,'_',stainfo(jsta).staname,'_',num2str(itime)]);
			filename = sprintf('data/%s_%d.mat',stainfo(ista).staname,itime);
			if exist(filename,'file')
				disp(['Cannot find ',filename]);
				datai = load(filename);
			else
				continue;
			end
			filename = sprintf('data/%s_%d.mat',stainfo(jsta).staname,itime);
			if exist(filename,'file')
				dataj = load(filename);
			else
				disp(['Cannot find ',filename]);
				continue;
			end
			xcor_sum=0;
			xcornum=0;
			for iseg = itime*TimeSegN:(itime+1)*TimeSegN-1
				if mod(iseg,100)==0
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
				end % Make sure station j is stable


				% Pre-process the data from two stations
				if iscor

					d1 = datai.segdata(iseg,:);
					d2 = datai.segdata(iseg,:);
					fftd1 = fft(d1);
					fftd2 = fft(d2);
					% Remove instrument response
					resp1 = cropfft(stainfo(ista).resp,3600);
					resp2 = cropfft(stainfo(jsta).resp,3600);
					fftd1 = fftd1.*resp1;
					fftd2 = fftd2.*resp2;
					% Whiten
					fftd1 = spectrumwhiten(fftd1);
					fftd2 = spectrumwhiten(fftd2);
					% inverse fft
					d1 = real(ifft(fftd1));
					d2 = real(ifft(fftd2));
					% Down sample to 1s
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
		end % end of station j
	end  % end of station i
end  % end of large time seg
