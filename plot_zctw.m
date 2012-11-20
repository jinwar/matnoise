load xspinfo_zc.mat

figure(123)
clf
hold on
for i = 1:length(xspinfo)
	if ~isempty(xspinfo(i).zctw)
		plot(twloc/2/pi,xspinfo(i).r./xspinfo(i).zctw);
	end
end
