function plot_ray_vel(ip)

load xspinfo_zc.mat
load seiscmap.mat
load coor.mat

global stainfo

r = 0.1;
rrange = [3 6];

for ixsp = 1:length(xspinfo)
	phv(ixsp) = xspinfo(ixsp).r./xspinfo(ixsp).tw(ip);
	rays(ixsp,1) = stainfo(xspinfo(ixsp).sta1).lat;
	rays(ixsp,2) = stainfo(xspinfo(ixsp).sta1).lon;
	rays(ixsp,3) = stainfo(xspinfo(ixsp).sta2).lat;
	rays(ixsp,4) = stainfo(xspinfo(ixsp).sta2).lon;
end
avgphv = nanmean(phv);

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
		v = xspinfo(ixsp).r./xspinfo(ixsp).tw(ip);
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

	figure(25)
	clf
	ax = worldmap(lalim, lolim);
	set(ax, 'Visible', 'off')
	for ixsp = 1:length(xspinfo)
		if xspinfo(ixsp).r < avgphv*rrange(1)*periods(ip) ||...
			xspinfo(ixsp).r > avgphv*rrange(2)*periods(ip)
			continue;
		end
		v = xspinfo(ixsp).r./xspinfo(ixsp).tw(ip);
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
end
