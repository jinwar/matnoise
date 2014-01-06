function dataout = runwin_norm(datain)
	N = 10;
	smamp = smooth(abs(datain),N);
	dataout = datain(:)./smamp(:);
return
