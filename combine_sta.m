% for CDpapua project, for station 13 and 32 turns to be the same station, all the xcor from these two stations should be stacked together.
%

clear;
load stainfo_BHZ.mat
load xspinfo.mat

% add station 32's data to station 13
ind = find([xspinfo(:).sta1] == 13 | [xspinfo(:).sta2]==13);

for ixsp = ind
	if xspinfo(ixsp).sta1 == 13
		sta2 = xspinfo(ixsp).sta2;
		data13 = load(xspinfo(ixsp).filename);
		cohere_sum = data13.cohere_sum;
		coherenum = data13.coherenum;
		ind32 = find([xspinfo(:).sta2] == sta2 & [xspinfo(:).sta1] == 32);
		if ~isempty(ind32)
			data32 = load(xspinfo(ind32).filename);
			cohere_sum = cohere_sum + data32.cohere_sum;
			coherenum  = coherenum + data32.coherenum;
			delete(xspinfo(ind32).filename);
			delfilename = xspinfo(ind32).filename;
		end
		ind32 = find([xspinfo(:).sta1] == sta2 & [xspinfo(:).sta2]== 32);
		if ~isempty(ind32)
			data32 = load(xspinfo(ind32).filename);
			cohere_sum = cohere_sum + conj(data32.cohere_sum);
			coherenum  = coherenum + data32.coherenum;
			delete(xspinfo(ind32).filename);
			delfilename = xspinfo(ind32).filename;
		end
	elseif xspinfo(ixsp).sta2 == 13
		sta1 = xspinfo(ixsp).sta1;
		data13 = load(xspinfo(ixsp).filename);
		cohere_sum = data13.cohere_sum;
		coherenum = data13.coherenum;
		ind32 = find([xspinfo(:).sta1] == sta1 & [xspinfo(:).sta2] == 32);
		if ~isempty(ind32)
			data32 = load(xspinfo(ind32).filename);
			cohere_sum = cohere_sum + data32.cohere_sum;
			coherenum  = coherenum + data32.coherenum;
			delete(xspinfo(ind32).filename);
			delfilename = xspinfo(ind32).filename;
		end
		ind32 = find([xspinfo(:).sta2] == sta1 & [xspinfo(:).sta1] == 32);
		if ~isempty(ind32)
			data32 = load(xspinfo(ind32).filename);
			cohere_sum = cohere_sum + conj(data32.cohere_sum);
			coherenum  = coherenum + data32.coherenum;
			delete(xspinfo(ind32).filename);
			delfilename = xspinfo(ind32).filename;
		end
	end
	save(xspinfo(ixsp).filename,'cohere_sum','coherenum')
	disp(['Saved: ',xspinfo(ixsp).filename,' Deleted: ',delfilename]);
end



