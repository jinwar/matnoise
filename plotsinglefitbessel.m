function plotsinglefitbessel(sta1,sta2)

global tN
global waxis
global twloc
global weight

load xspinfo
load stainfo_BHZ

frange = [0.04 0.15];
tN = 20;
refc1 = 3.6;
refc2 = 2.5;
Isfigure = 0;

refc = refc1 + (refc2-refc1)/(tN)*[1:tN];

twloc = frange(1):(frange(2)-frange(1))/(tN-1):frange(2);
twloc = twloc*2*pi;
waxis = (frange(1):1/3600:frange(2))*2*pi;
goodnum=0;
weight = waxis;
weight(:) = 1;
ind = find([xspinfo(:).sta1] == sta1 & [xspinfo(:).sta2] == sta2);
if isempty(ind)
	disp('Can not find station pair');
	return
end

besselerr(xspinfo(ind).tw,xspinfo(ind).xsp,1);
figure(1)
title(num2str(xspinfo(ind).sumerr));
subplot(3,1,1)
title([stainfo(sta1).staname,'--',stainfo(sta2).staname]);
figure(2)
% clf
hold on
faxis = 2*pi./twloc;
phv = xspinfo(ind).r./xspinfo(ind).tw;
plot(faxis,phv,'b-o');
text(faxis(1),phv(1),...
    [stainfo(sta1).staname,'--',stainfo(sta2).staname]);
        
end
