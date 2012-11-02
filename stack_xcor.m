% Matlab Script to read in the cross-correlation mat data files under the folder ./xcor/
% and stack them based on station pairs.
% Written by Ge Jin, oct 2012

clear 
load stainfo_BHZ.mat
Isfigure=1;
for ista = 1:length(stainfo)
% for ista = 1:1
	for jsta = ista+1:length(stainfo)
%     for jsta = 4:4
		xcor_sum = 0;
		xcornum = 0;
		disp(['stacking ',stainfo(ista).staname,' and ',stainfo(jsta).staname]);
		for itime = 1:11
			filename = sprintf('xcor/%s_%s_%d.mat',stainfo(ista).staname,stainfo(jsta).staname,itime);
			if ~exist(filename,'file')
				disp(['Cannot find ', filename]);
				continue;
			end
			data = load(filename);
			xcor_sum = xcor_sum + data.xcor_sum;
			xcornum = xcornum + data.xcornum;
		end % loop of itime
		if xcornum > 1
			xcor_avg = xcor_sum./xcornum;
			filename = sprintf('xcor/%s_%s.mat',stainfo(ista).staname,stainfo(jsta).staname);
			save(filename,'xcor_avg','xcornum');
		else
			disp(['No xcor for ',stainfo(ista).staname,' and ',stainfo(jsta).staname]);
		end
	end % loop of jsta
end % loop of ista
