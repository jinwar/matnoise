function plot_avg_tomo(ip)

load eikonal_avg.mat
load seiscmap
r=0.2;

figure(40)
clf
lalim = [min(xnode) max(xnode)];
lolim = [min(ynode) max(ynode)];
[xi yi] = ndgrid(xnode,ynode);
ax = worldmap(lalim, lolim);
set(ax, 'Visible', 'off')
surfacem(xi,yi,avgtomo(ip).GV);
drawpng
title(['Periods: ',num2str(ip)],'fontsize',15)
avgv = nanmean(avgtomo(ip).GV(:));
caxis([avgv*(1-r) avgv*(1+r)])
colorbar
colormap(seiscmap)

figure(41)
clf
lalim = [min(xnode) max(xnode)];
lolim = [min(ynode) max(ynode)];
[xi yi] = ndgrid(xnode,ynode);
ax = worldmap(lalim, lolim);
set(ax, 'Visible', 'off')
surfacem(xi,yi,avgtomo(ip).raydense);
drawpng
title(['Ray Dense Periods: ',num2str(ip)],'fontsize',15)
colorbar

end
