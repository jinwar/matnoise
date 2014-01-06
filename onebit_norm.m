function outdata = onebit_norm(indata)
	outdata = detrend(indata);
	outdata(find(indata > 0)) = 1;
	outdata(find(indata < 0)) = -1;
end
