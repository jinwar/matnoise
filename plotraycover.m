clear

load stainfo_BHZ
load xspinfo
lalim = [-10.8 -8.2];
lolim = [148.8 151.5];
figure(5)
clf
hold on
ax = worldmap(lalim, lolim);
set(ax, 'Visible', 'off')

goodnum=0;
for ind =1:length(xspinfo)
    if xspinfo(ind).sumerr <1
        dx = diff(xspinfo(ind).tw.*twloc);
        badind = find(dx<0);
        if length(badind)==0
            goodnum = goodnum+1;
            lat(1) = stainfo(xspinfo(ind).sta1).lat;
            lat(2) = stainfo(xspinfo(ind).sta2).lat;
            lon(1) = stainfo(xspinfo(ind).sta1).lon;
            lon(2) = stainfo(xspinfo(ind).sta2).lon;
			plotm(lat,lon,'b');
        end
    end
end

drawpng

for ista = 1:length(stainfo)
	plotm(stainfo(ista).lat,stainfo(ista).lon,'rv')
	textm(stainfo(ista).lat+0.05,stainfo(ista).lon+0.05,stainfo(ista).staname,'color','r');
end
