global tN
global waxis
global twloc
global weight

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

for ind =1:length(xspinfo)
    if xspinfo(ind).sumerr <1
        dx = diff(xspinfo(ind).tw.*twloc);
        badind = find(dx<0);
        if length(badind)==0
            besselerr(xspinfo(ind).tw,xspinfo(ind).xsp,1);
            figure(1)
            title(num2str(xspinfo(ind).sumerr));
            figure(2)
            hold on
            plot(twloc/2/pi,xspinfo(ind).r./xspinfo(ind).tw);
            goodnum = goodnum+1;
% 			pause(1)
%             pause
        end
        
    end
end
goodnum
