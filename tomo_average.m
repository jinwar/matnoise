clear

load event_tomo.mat
load seiscmap.mat

r = 0.2;

[m n] = size(event_tomo(1).GV);
waterlevel = [2 4];
for ip = 1:size(event_tomo,2)
	sumGV = zeros(size(event_tomo(1).GV));
	sumweight = zeros(size(event_tomo(1).GV));
	for ie = 1:size(event_tomo,1)
		if size(event_tomo(ie,ip).GV,1)~=m
			continue;
		end
		for ix = 1:m
			for iy = 1:n
				ind = find(event_tomo(ie,ip).GV>waterlevel(2));
				event_tomo(ie,ip).GV(ind) = waterlevel(2);
				ind = find(event_tomo(ie,ip).GV<waterlevel(1));
				event_tomo(ie,ip).GV(ind) = waterlevel(1);
				if ~isnan(event_tomo(ie,ip).GV(ix,iy))
					sumGV(ix,iy) = sumGV(ix,iy) + event_tomo(ie,ip).GV(ix,iy)*event_tomo(ie,ip).raydense(ix,iy);
					sumweight(ix,iy) = sumweight(ix,iy) + event_tomo(ie,ip).raydense(ix,iy);
%					sumGV(ix,iy) = sumGV(ix,iy) + event_tomo(ie,ip).GV(ix,iy);
%					sumweight(ix,iy) = sumweight(ix,iy) + 1;
				end
			end
		end
	end
	avgGV = sumGV./sumweight;
	average_tomo(ip).GV = avgGV;
	avgtomo(ip).GV = avgGV;
	avgtomo(ip).raydense = sumweight;
	avgtomo(ip).avgV = nanmean(avgGV(:));
end

% Calculate the variance
for ip=1:size(event_tomo,2)
	sumvar = zeros(size(event_tomo(1).GV));
	for ie = 1:size(event_tomo,1)
		if size(event_tomo(ie,ip).GV,1)~=m
			continue;
		end
		for ix = 1:m
			for iy = 1:n
				if ~isnan(event_tomo(ie,ip).GV(ix,iy))
					sumvar(ix,iy) = sumvar(ix,iy) + ...
						(event_tomo(ie,ip).GV(ix,iy)-avgtomo(ip).GV(ix,iy)).^2*event_tomo(ie,ip).raydense(ix,iy);
				end
			end
		end
	end
	for ix = 1:m
		for iy = 1:n
			if isnan(avgtomo(ip).GV(ix,iy))
				sumvar(ix,iy) = NaN;
			end
		end
	end
	avgtomo(ip).GVvar = sqrt(sumvar./sumweight);
end

save('eikonal_avg.mat','avgtomo','xnode','ynode','periods');

figure(15)
clf
lalim = [min(xnode) max(xnode)];
lolim = [min(ynode) max(ynode)];
[xi yi] = ndgrid(xnode,ynode);
for ip = 1:size(event_tomo,2)
    subplot(4,5,ip)
    ax = worldmap(lalim, lolim);
    set(ax, 'Visible', 'off')
    surfacem(xi,yi,avgtomo(ip).GV);
    drawpng
    title(['Periods: ',num2str(periods(ip))],'fontsize',15)
	avgv = nanmean(avgtomo(ip).GV(:));
    caxis([avgv*(1-r) avgv*(1+r)])
    colorbar
	colormap(seiscmap)
end
figure(16)
clf
lalim = [min(xnode) max(xnode)];
lolim = [min(ynode) max(ynode)];
[xi yi] = ndgrid(xnode,ynode);
for ip = 1:size(event_tomo,2)
    subplot(4,5,ip)
    ax = worldmap(lalim, lolim);
    set(ax, 'Visible', 'off')
    surfacem(xi,yi,avgtomo(ip).GVvar);
    drawpng
    title(['Periods: ',num2str(periods(ip))],'fontsize',15)
	colorbar
	caxis([0 1])
end
figure(17)
clf
lalim = [min(xnode) max(xnode)];
lolim = [min(ynode) max(ynode)];
[xi yi] = ndgrid(xnode,ynode);
for ip = 1:size(event_tomo,2)
    subplot(4,5,ip)
    ax = worldmap(lalim, lolim);
    set(ax, 'Visible', 'off')
    surfacem(xi,yi,avgtomo(ip).raydense);
    drawpng
    title(['Periods: ',num2str(periods(ip))],'fontsize',15)
    colorbar
end
