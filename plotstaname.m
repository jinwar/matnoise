clear

load stainfo_BHZ
load xspinfo
lalim = [-10.8 -8.2];
lolim = [148.8 151.5];
figure(5)
clf
hold on
ax = worldmap(lalim, lolim);
load pngcoastline
geoshow([S.Lat], [S.Lon], 'Color', [0.5 0.5 0.5],'linewidth',1)
set(ax, 'Visible', 'off')

for ista = 1:length(stainfo)
	plotm(stainfo(ista).lat,stainfo(ista).lon,'rv')
	textm(stainfo(ista).lat+0.05,stainfo(ista).lon,stainfo(ista).staname,'color','k');
	textm(stainfo(ista).lat,stainfo(ista).lon+0.03,num2str(ista),'color','k');
end
