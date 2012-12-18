load raytomo.mat
load eikonal_avg.mat
load coor.mat
r = 0.2;
ip=8;
avgphv = nanmean(raytomo(ip).GV(:));

figure(51)
clf
ax = worldmap(lalim, lolim);
set(ax, 'Visible', 'off');
surfacem(xi,yi,avgtomo(ip).GV);
load seiscmap
colormap(seiscmap)
drawpng
colorbar
caxis([avgphv*(1-r) avgphv*(1+r)]);
title('Eikonal');
figure(52)
clf
ax = worldmap(lalim, lolim);
set(ax, 'Visible', 'off');
surfacem(xi,yi,raytomo(ip).GV);
load seiscmap
colormap(seiscmap)
drawpng
colorbar
caxis([avgphv*(1-r) avgphv*(1+r)]);
title('Ray');
figure(53)
clf
ax = worldmap(lalim, lolim);
set(ax, 'Visible', 'off');
surfacem(xi,yi,avgtomo(ip).GV-raytomo(ip).GV);
load seiscmap
colormap(seiscmap)
drawpng
colorbar
disp(['Error to Ray tomo: ' ,num2str(nanmean(abs(avgtomo(ip).GV(:)-raytomo(ip).GV(:))))]);