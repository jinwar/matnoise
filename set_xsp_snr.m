clear

load stainfo_BHZ.mat
load xspinfo.mat

faxis = [0:1799]/3600;
signalind = find(faxis > 0.04 & faxis < 0.15);
noiseind = find(faxis > 0.4 );

for ixsp = 1:length(xspinfo)
	sta1 = xspinfo(ixsp).sta1;
	sta2 = xspinfo(ixsp).sta2;
	filename = sprintf('xcorf/%s_%s_f.mat',stainfo(sta1).staname,stainfo(sta2).staname);
	data1 = load(filename);
	xcorf1 = data1.cohere_sum./data1.coherenum;
	signal = mean(abs((xcorf1(signalind))).^1);
	noise = mean(abs((xcorf1(noiseind))).^1);
	snr = signal/noise;
	xspinfo(ixsp).filename = filename;
	xspinfo(ixsp).snr = snr;
end

save('xspinfo.mat','twloc','waxis','xspinfo');
