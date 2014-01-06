% This script is used to calculate the frequency domain ambient noise using Aki's method

clear

% Parameter settings
StableSeg = 1;   % How many hour the instrument stable after being powered
TimeSegN = 1000;
% system
Isoverwrite = 0;
Isfigure = 0;
isoutput = 1;
% Normalization
is_time_norm = 1;
is_spectrum_whiten = 0;
is_remove_events = 0;

load stainfo_BHZ.mat


%
 for ista = 1:length(stainfo)-1
     for jsta = ista+1:length(stainfo)
%for ista = 1
%    for jsta = 2:3
        % Check existing datafile for recovering from last run
        filename = ['xcorf/',stainfo(ista).staname,'_',stainfo(jsta).staname,'_f.mat'];
        if exist(filename,'file') && ~Isoverwrite
            disp(['Exist ',filename,' Skip!']);
            continue;
        end
        % start to loop through the times
        cohere_sum=0;
        coherenum=0;
        for itime = 0:11
            
            disp([stainfo(ista).staname,'_',stainfo(jsta).staname,'_',num2str(itime)]);
            filename = sprintf('data/%s_%d.mat',stainfo(ista).staname,itime);
            if exist(filename,'file')
                datai = load(filename);
            else
                disp(['Cannot find ',filename]);
                continue;
            end
            filename = sprintf('data/%s_%d.mat',stainfo(jsta).staname,itime);
            if exist(filename,'file')
                dataj = load(filename);
            else
                disp(['Cannot find ',filename]);
                continue;
            end
            for iseg = itime*TimeSegN:(itime+1)*TimeSegN-1
                if iseg == 0
                    continue;
                end
                % Test whether we want to make the xcor
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
                if iseg > size(datai.segdata,1) || iseg > size(dataj.segdata,1)
                    iscor = 0;
                end
                
                
                % Pre-process the data from two stations
                if iscor
                    d1 = datai.segdata(iseg,:);
                    d2 = dataj.segdata(iseg,:);
                    d1 = detrend(d1);
                    d2 = detrend(d2);
					if is_time_norm
						d1 = runwin_norm(d1);
						d2 = runwin_norm(d2);
					end
                    fftd1 = fft(d1);
                    fftd2 = fft(d2);
                    % Remove instrument response
                    resp1 = cropfft(stainfo(ista).resp,3600);
                    resp2 = cropfft(stainfo(jsta).resp,3600);
                    fftd1 = fftd1(:).*resp1(:);
                    fftd2 = fftd2(:).*resp2(:);
                    % Whiten
                    if is_spectrum_whiten
                        fftd1 = spectrumwhiten(fftd1);
                        fftd2 = spectrumwhiten(fftd2);
                    end
                    % inverse fft
                    d1 = real(ifft(fftd1));
                    d2 = real(ifft(fftd2));
                    if isnan(sum(d2)) || isnan(sum(d1))
                        iscor=0;
                    end
                    % Down sample to 1s
                    if length(d1)~=length(d2)
                        disp(['Timeseg: ',num2str(iseg), ...
                            'Down Sample Length error' ]);
                        iscor=0;
                    end
                end
                
                % Do Xcor!
                if iscor
                    this_cohere = fftd1.*conj(fftd2);
                    this_cohere=this_cohere./abs(fftd1)./abs(fftd2);
                    cohere_sum = cohere_sum + this_cohere;
                    coherenum = coherenum+1;
                end
            end % Loop of time segments
        end  % end of large time seg
        if coherenum > 1
            if Isfigure
                figure
                clf
				set(gcf,'position',[400 400 600 300]);
                dt = 1;
                T = length(cohere_sum);
                faxis = [0:1/T:1/dt/2,-1/dt/2+1/T:1/T:-1/T];
                ind = find(faxis>0);
                plot(faxis(ind),smooth(real(cohere_sum(ind)/coherenum),10));
                dist = distance(stainfo(ista).lat,stainfo(ista).lon,stainfo(jsta).lat,stainfo(jsta).lon);
                dist = deg2km(dist);
                title(sprintf('%s %s coherency,station distance: %f km',stainfo(ista).staname,stainfo(jsta).staname,dist));
                xlim([0.04 0.16])
				drawnow
                %				pause
            end
            if isoutput
                save(sprintf('xcorf/%s_%s_f.mat',stainfo(ista).staname,stainfo(jsta).staname),'cohere_sum','coherenum');
            end
        end % end of if coherenum
    end % end of station j
end  % end of station i

