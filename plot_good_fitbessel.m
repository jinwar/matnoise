clear;

load xspinfo.mat

goodind = find([xspinfo(:).r]<80 & [xspinfo(:).sumerr]<0.5 & ...
    [xspinfo(:).coherenum] > 4000);

% for ixsp = goodind
% 	ixsp
% 	sta1 = xspinfo(ixsp).sta1;
% 	sta2 = xspinfo(ixsp).sta2;
% 	plotsinglefitbessel(sta1,sta2); 
% 	pause
% end
% 
frange = [0.04 0.15];
waxis = (frange(1):1/3600:frange(2))*2*pi;

sta1 = 1; sta2 = 2;
ixsp1 = find([xspinfo(:).sta1]==sta1 & [xspinfo(:).sta2]==sta2)
sta1 = 1; sta2 = 3;
ixsp2 = find([xspinfo(:).sta1]==sta1 & [xspinfo(:).sta2]==sta2)
%ixsp2 = 728;
fontsize =15;

sta1 = xspinfo(ixsp1).sta1;
sta2 = xspinfo(ixsp1).sta2;
tw = xspinfo(ixsp1).tw;
xsp = xspinfo(ixsp1).xsp;
tw1 = interp1(twloc,tw,waxis,'linear'); 
x1 = waxis.*tw1;
be = besselj(0,x1);
be = be./mean(abs(be)).*mean([abs(xsp)]);
figure(75)
clf
subplot(2,1,1)
hold on
plot(waxis/2/pi,xsp,'b','linewidth',2);
plot(waxis/2/pi,be,'b--','linewidth',2);
xlabel('Frequency (Hz)','fontsize',fontsize);
ylabel('Amplitude','fontsize',fontsize);
sta1 = xspinfo(ixsp2).sta1;
sta2 = xspinfo(ixsp2).sta2;
tw = xspinfo(ixsp2).tw;
xsp = xspinfo(ixsp2).xsp;
tw1 = interp1(twloc,tw,waxis,'linear'); 
x1 = waxis.*tw1;
be = besselj(0,x1);
be = be./mean(abs(be)).*mean([abs(xsp)]);
subplot(2,1,2)
hold on
plot(waxis/2/pi,xsp,'r','linewidth',2);
plot(waxis/2/pi,be,'r--','linewidth',2);
xlabel('Frequency (Hz)','fontsize',fontsize);
ylabel('Amplitude','fontsize',fontsize);
set(gcf,'position',[100 300 700 400])
%export_fig('pics/besselfit/besselfitwaveform','-transparent','-m2');

figure(76)
clf
hold on
plot(2*pi./twloc,xspinfo(ixsp1).r./xspinfo(ixsp1).tw,'bo-','linewidth',2,'markersize',10)
plot(2*pi./twloc,xspinfo(ixsp2).r./xspinfo(ixsp2).tw,'ro-','linewidth',2,'markersize',10)
xlabel('Period (s)','fontsize',fontsize);
ylabel('Phase Velocity (km/s)','fontsize',fontsize);
set(gca,'fontsize',fontsize);
set(gcf,'position',[800 300 300 400])
%export_fig('pics/besselfit/besselfitphv','-transparent','-m2');

