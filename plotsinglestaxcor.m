function plotsinglestaxcor(ista)
% Matlab function to plot all the xcor with one station again station stadistance
% Written by Ge Jin, jinwar@gmail.com
%

% ista = 1
coperiod = [5 20];

f1 = 1/coperiod(2);
f2 = 1/coperiod(1);
[b a] = butter(2,[f1 f2]);

load stainfo_BHZ.mat

stapairn = 0;
for jsta = 1:length(stainfo)
	if jsta < ista
		filename = sprintf('xcor/%s_%s.mat',stainfo(jsta).staname,stainfo(ista).staname);
	else
		filename = sprintf('xcor/%s_%s.mat',stainfo(ista).staname,stainfo(jsta).staname);
	end
	if exist(filename,'file')
		disp(['Read in ',filename])
		data =  load(filename);
		stapairn = stapairn+1;
		stapairxcor(stapairn,:) = filtfilt(b,a,data.xcor_avg);
		stapairxcor(stapairn,:) = stapairxcor(stapairn,:)/max(abs(stapairxcor(stapairn,:)));
		stadist(stapairn) = distance(stainfo(ista).lat,stainfo(ista).lon,stainfo(jsta).lat,stainfo(jsta).lon);
        staid(stapairn) = jsta;
	end
end
if stapairn ==0
    return
end
stadist = deg2km(stadist);
N = length(data.xcor_avg);
time = [-(N-1)/2:(N-1)/2];

amp = 1e1;
figure(1)
clf
hold on
for istap = 1:stapairn
	plot(time,stapairxcor(istap,:)*amp+stadist(istap),'k')
    text(150,stadist(istap),stainfo(staid(istap)).staname);
end
xlim([-200 200])
title(stainfo(ista).staname);


		

