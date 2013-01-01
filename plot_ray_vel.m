function plot_ray_vel(ip)

load xspinfo.mat
load seiscmap.mat
load coor.mat
load stainfo_BHZ.mat
load raytomo.mat


r = 0.2;
rrange = [1 6];
if ~exist('ip','var')
    ip=16;
end

errlevel = 1;
snrtol = 1.1;
% Initial the xsp structure
for ixsp = 1:length(xspinfo)
	xspinfo(ixsp).isgood = 0;
	if xspinfo(ixsp).sumerr < errlevel ...
            && xspinfo(ixsp).snr > snrtol && xspinfo(ixsp).coherenum > 2000
        xspinfo(ixsp).isgood = 1;
	end
end % end of loop ixsp

for ixsp = 1:length(xspinfo)
    phv(ixsp) = xspinfo(ixsp).r./xspinfo(ixsp).tw(ip);
    rays(ixsp,1) = stainfo(xspinfo(ixsp).sta1).lat;
    rays(ixsp,2) = stainfo(xspinfo(ixsp).sta1).lon;
    rays(ixsp,3) = stainfo(xspinfo(ixsp).sta2).lat;
    rays(ixsp,4) = stainfo(xspinfo(ixsp).sta2).lon;
end
avgphv = nanmean(raytomo(ip).GV(:));

periods = 2*pi./twloc;
vrange = avgphv.*[1-r,1+r];
vx = linspace(vrange(1),vrange(2),size(seiscmap,1));

figure(24)
clf
ax = worldmap(lalim, lolim);
set(ax, 'Visible', 'off')
for ixsp = 1:length(xspinfo)
    if xspinfo(ixsp).r < avgphv*rrange(1)*periods(ip) ||...
            xspinfo(ixsp).r > avgphv*rrange(2)*periods(ip)
        continue;
    end
    if xspinfo(ixsp).isgood==0
        continue;
    end
    v = phv(ixsp);
    if v<vrange(1)
        v = vrange(1);
    elseif v>vrange(2)
        v = vrange(2);
    end
    linecolor = interp1(vx,seiscmap,v);
    lat(1) = stainfo(xspinfo(ixsp).sta1).lat;
    lon(1) = stainfo(xspinfo(ixsp).sta1).lon;
    lat(2) = stainfo(xspinfo(ixsp).sta2).lat;
    lon(2) = stainfo(xspinfo(ixsp).sta2).lon;
    plotm(lat,lon,'color',linecolor);
end
for ista = 1:length(stainfo)
    plotm(stainfo(ista).lat,stainfo(ista).lon,'kv','markersize',10,'MarkerFaceColor','k')
end
drawpng
colormap(seiscmap);
colorbar
caxis(vrange);
title(['Periods: ',num2str(periods(ip))],'fontsize',15)


figure(25)
clf
ax = worldmap(lalim, lolim);
set(ax, 'Visible', 'off')
for ixsp = 1:length(xspinfo)
    if xspinfo(ixsp).r < avgphv*rrange(1)*periods(ip) ||...
            xspinfo(ixsp).r > avgphv*rrange(2)*periods(ip)
        continue;
    end
    v = phv(ixsp);
    if v<vrange(1)
        v = vrange(1);
    elseif v>vrange(2)
        v = vrange(2);
    end
    linecolor = interp1(vx,seiscmap,v);
    lat(1) = stainfo(xspinfo(ixsp).sta1).lat;
    lon(1) = stainfo(xspinfo(ixsp).sta1).lon;
    lat(2) = stainfo(xspinfo(ixsp).sta2).lat;
    lon(2) = stainfo(xspinfo(ixsp).sta2).lon;
    [midlat midlon] = gcwaypts(lat(1),lon(1),lat(2),lon(2),2);
    plotm(midlat(2),midlon(2),'.','color',linecolor,'markersize',20);
end
for ista = 1:length(stainfo)
    plotm(stainfo(ista).lat,stainfo(ista).lon,'kv','markersize',10,'MarkerFaceColor','k')
end
drawpng
colormap(seiscmap);
colorbar
caxis(vrange);
title(['Periods: ',num2str(periods(ip))],'fontsize',15)

figure(26)
clf
ax = worldmap(lalim, lolim);
set(ax, 'Visible', 'off')
surfacem(xi,yi,raytomo(ip).GV);
caxis([avgphv*(1-r) avgphv*(1+r)]);
drawpng
colormap(seiscmap);
colorbar
% 	plot_avg_tomo(ip);
end
