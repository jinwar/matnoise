load raytomo;
load seiscmap
r=0.2;


for ip = 1:20
    figure(40)
    clf
    lalim = [min(xnode) max(xnode)];
    lolim = [min(ynode) max(ynode)];
    [xi yi] = ndgrid(xnode,ynode);
    ax = worldmap(lalim, lolim);
    set(ax, 'Visible', 'off')
    axesm(gcm,'fontsize',10);
    h1=surfacem(xi,yi,raytomo(ip).GV);
    % set(h1,'facecolor','interp');
    drawpng
    avgv = nanmean(raytomo(ip).GV(:));
    caxis([avgv*(1-r) avgv*(1+r)])
    colorbar
    colormap(seiscmap)
    filename = sprintf('pics/raytomo/phasev_%2d',ip);
    filename(find(filename==' '))='0';
%     export_fig(filename,'-transparent','-png','-m2')
end
