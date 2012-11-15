function plotcsray(event,stainfo,ie)

	lalim = [-10.8 -8.2];
	lolim = [148.8 151.5];
	figure(5)
	clf
	for ip=1:size(event,2)
		evla = event(ie,ip).evla;
		evlo = event(ie,ip).evlo;
		ray = event(ie,ip).ray;
		subplot(4,5,ip)
		ax = worldmap(lalim, lolim);
		set(ax, 'Visible', 'off')
		for ind =1:size(ray,1)
			lat(1) = ray(ind,1);
			lat(2) = ray(ind,3);
			lon(1) = ray(ind,2);
			lon(2) = ray(ind,4);
			plotm(lat,lon,'b');
		end
		for ista = 1:length(stainfo)
			plotm(stainfo(ista).lat,stainfo(ista).lon,'rv')
	%		textm(stainfo(ista).lat+0.05,stainfo(ista).lon+0.05,stainfo(ista).staname,'color','r');
		end
		plotm(evla,evlo,'ko','markersize',15)
		drawpng
	end

end
