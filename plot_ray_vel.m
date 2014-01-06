function plot_ray_vel(ip)

load xspinfo.mat
load seiscmap.mat
load coor.mat
load stainfo_BHZ.mat
load raytomo.mat


setup_parameters;

if ~exist('ip','var')
    ip=16;
end

ip

% read in bad station list, if existed
if exist('badsta.lst')
	badstnms = textread('badsta.lst','%s');
	badstaids = find(ismember({stainfo.staname},badstnms));
	disp('Found Bad stations:')
	disp(badstnms)
end

% Initial the xsp structure
for ixsp = 1:length(xspinfo)
	xspinfo(ixsp).isgood = 0;
	if xspinfo(ixsp).sumerr < errlevel ...
            && xspinfo(ixsp).snr > snrtol && xspinfo(ixsp).coherenum > 2000
        xspinfo(ixsp).isgood = 1;
	end
	if sum(ismember([xspinfo(ixsp).sta1 xspinfo(ixsp).sta2],badstaids)) > 0
		xspinfo(ixsp).isgood = 0;
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
    if xspinfo(ixsp).r < refv*distrange(1)*periods(ip) ||...
            xspinfo(ixsp).r > refv*distrange(2)*periods(ip)
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
    if xspinfo(ixsp).r < refv*distrange(1)*periods(ip) ||...
            xspinfo(ixsp).r > refv*distrange(2)*periods(ip)
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
for ista = 1:length(stainfo)
    plotm(stainfo(ista).lat,stainfo(ista).lon,'kv','markersize',10,'MarkerFaceColor','k')
end
% 	plot_avg_tomo(ip);

figure(27)
clf
ax = worldmap(lalim, lolim);
set(ax, 'Visible', 'off')
err_range = [-3 3];
errx = linspace(err_range(1),err_range(2),size(seiscmap,1));
for ixsp = 1:size(raytomo(ip).rays,1)
%	if raytomo(ip).w(ixsp) == 0
%		continue;
%	end
    err = raytomo(ip).err(ixsp);
    if err<err_range(1)
        err = err_range(1);
    elseif err>err_range(2)
        err = err_range(2);
    end
    linecolor = interp1(errx,seiscmap,err);
    lat(1) = raytomo(ip).rays(ixsp,1);
    lon(1) = raytomo(ip).rays(ixsp,2);
    lat(2) = raytomo(ip).rays(ixsp,3);
    lon(2) = raytomo(ip).rays(ixsp,4);
    plotm(lat,lon,'color',linecolor);
end
for ista = 1:length(stainfo)
    plotm(stainfo(ista).lat,stainfo(ista).lon,'kv','markersize',10,'MarkerFaceColor','k')
end
drawpng
colormap(seiscmap);
colorbar
caxis(err_range);
title(['Periods: ',num2str(periods(ip))],'fontsize',15)

end
