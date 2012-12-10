function plot_avg_tomo(ip)

load eikonal_avg.mat
load seiscmap
load xspinfo.mat
r=0.2;

periods = 2*pi./twloc;

figure(40)
clf
lalim = [min(xnode) max(xnode)];
lolim = [min(ynode) max(ynode)];
[xi yi] = ndgrid(xnode,ynode);
ax = worldmap(lalim, lolim);
set(ax, 'Visible', 'off')
h1=surfacem(xi,yi,avgtomo(ip).GV);
% set(h1,'facecolor','interp');
drawpng
title(['Periods: ',num2str(periods(ip))],'fontsize',15)
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
title(['Ray Dense Periods: ',num2str(periods(ip))],'fontsize',15)
colorbar

figure(42)
clf
lalim = [min(xnode) max(xnode)];
lolim = [min(ynode) max(ynode)];
[xi yi] = ndgrid(xnode,ynode);
ax = worldmap(lalim, lolim);
set(ax, 'Visible', 'off')
surfacem(xi,yi,avgtomo(ip).GVvar);
drawpng
title(['Variance Periods: ',num2str(periods(ip))],'fontsize',15)
caxis([0 1])
colorbar

end
