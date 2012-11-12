
lalim = [-10.8 -8.2];
lolim = [148.8 151.5];
figure(5)
clf
hold on
ax = worldmap(lalim, lolim);
set(ax, 'Visible', 'off')

goodnum=0;
for ind =1:size(ray,1)
	lat(1) = ray(ind,1);
	lat(2) = ray(ind,3);
	lon(1) = ray(ind,2);
	lon(2) = ray(ind,4);
	plotm(lat,lon,'b');
end

drawpng

for ista = 1:length(stainfo)
	plotm(stainfo(ista).lat,stainfo(ista).lon,'rv')
	textm(stainfo(ista).lat+0.05,stainfo(ista).lon+0.05,stainfo(ista).staname,'color','r');
end

plotm(evla,evlo,'ko','markersize',15)