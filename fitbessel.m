% Program to fit the xcorf and get the travel time for each frequency for two station pairs

clear

global tN
global waxis
global twloc
global weight

load stainfo_BHZ.mat
load refphasev.mat
frange = [0.04 0.15];
tN = 20;
Isfigure =0;
isoutput = 1;
nearstadist = 5;

refc1 = 3.6;
refc2 = 2.7;
refc = refc1 + (refc2-refc1)/(tN)*[1:tN];
% refc = refv;

twloc = frange(1):(frange(2)-frange(1))/(tN-1):frange(2);
twloc = twloc*2*pi;
waxis = (frange(1):1/3600:frange(2))*2*pi;
faxis = [0:1800-1]*1/3600;
signalind = find(faxis > frange(1) & faxis < frange(2));
noiseind = find(faxis > 0.4 );

stapairn = 0;
sumerr = 0;
for stai = 1:length(stainfo)
% for stai = 1
    sta1 = stai;
       for staj = sta1+1:length(stainfo)
%     for staj = 14
        sta2 = staj;
        r1 = distance(stainfo(sta1).lat,stainfo(sta1).lon,stainfo(sta2).lat,stainfo(sta2).lon);
        r1 = deg2km(r1);
        if r1 < nearstadist
            continue;
        end
        filename = sprintf('xcorf/%s_%s_f.mat',stainfo(sta1).staname,stainfo(sta2).staname);
        if ~exist(filename,'file')
            disp(['not exist ',filename])
            continue;
        end
        data1 = load(filename);
        xcorf1 = data1.cohere_sum./data1.coherenum;
        
        N = 1000;
        xcorf1 = real(xcorf1(1:N));
        xcorf1(1) = 0;
        %xcorf1 = smooth(xcorf1,10);
        %xcorf2 = smooth(xcorf2,10);
        faxis = [0:N-1]*1/3600;
        
        
        xsp1 = interp1(faxis*2*pi,xcorf1,waxis);
        
        xsp1 = smooth(xsp1,10);
        
        tw1 = ones(1,tN)*r1./refc;
        
        weight  = 1./waxis;
        tw2 = lsqnonlin(@(x) besselerr(x,[xsp1]),[tw1],[tw1]*0.8,[tw1]*1.2);
        weight(:) = 1;
        tw = lsqnonlin(@(x) besselerr(x,[xsp1]),[tw2],[tw2]*0.8,[tw2]*1.2);
        tw=tw2;
        
        stapairn = stapairn+1;
        xspinfo(stapairn).sta1 = sta1;
        xspinfo(stapairn).sta2 = sta2;
        xspinfo(stapairn).r = r1;
        xspinfo(stapairn).tw = tw;
        xspinfo(stapairn).xsp = xsp1;
        xspinfo(stapairn).coherenum = data1.coherenum;
        err = besselerr(tw,xsp1);
        err = err(1:length(waxis));
        xspinfo(stapairn).sumerr = sum(err.^2)./sum((xsp1./weight(:)).^2);
        xspinfo(stapairn).err = err./weight(:);
        sumerr = sumerr + (xspinfo(stapairn).err);
        data(stapairn,:) = r1./tw;
        
		xcorf1 = data1.cohere_sum./data1.coherenum;
		signal = mean(abs((xcorf1(signalind))).^1);
		noise = mean(abs((xcorf1(noiseind))).^1);
		snr = signal/noise;
		xspinfo(stapairn).filename = filename;
		xspinfo(stapairn).snr = snr;

        disp([filename,' fitted'])
        if Isfigure
            besselerr(tw,xsp1,3);
            figure(2)
            clf
            hold on
            plot(twloc/2/pi,r1./tw1,'ro-');
            plot(twloc/2/pi,r1./tw,'ko-');
            %             pause
        end
        %     pause
    end %end of station j
end  %end of station i
stapairn
%soundsc(rand(2000,1),1000,8)
if isoutput
    save('xspinfo.mat','xspinfo','twloc','waxis');
end
