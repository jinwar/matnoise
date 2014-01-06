clear;

load stainfo_BHZ.mat
load coor.mat

stlas = [stainfo.lat];
stlos = [stainfo.lon];
stnms = {stainfo.staname};

for ista = 1:length(stlas)
	disp(sprintf('%d %f %f %s',ista,stlas(ista),stlos(ista),char(stnms(ista))));
	stids(ista) = {num2str(ista)};
end

figure(78)
clf
ax = worldmap(lalim, lolim);
set(ax, 'Visible', 'off')
axesm(gcm,'fontsize',15);
load pngcoastline
geoshow([S.Lat], [S.Lon], 'Color', 'blue','linewidth',1)
plotm(stlas,stlos,'rv');
textm(stlas,stlos+0.05,char(stnms));
textm(stlas-0.05,stlos,char(stids));

