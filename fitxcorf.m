% Program to fit the xcorf and get the travel time for each frequency for two station pairs

clear

global tN
global waxis
global twloc
global bzeros
global xcommon

load stainfo_BHZ.mat
frange = [0.04 0.15];
tN = 10;
refc1 = 3.6;
refc2 = 3.0;
refc = refc1 + (refc2-refc1)/(tN)*[1:tN];

twloc = frange(1):(frange(2)-frange(1))/(tN-1):frange(2);
twloc = twloc*2*pi;
bzeros = besselzero(0,100,1);

sta1 = 1;
sta2 = 14;
sta3 = 1;
sta4 = 25;

r1 = distance(stainfo(sta1).lat,stainfo(sta1).lon,stainfo(sta2).lat,stainfo(sta2).lon);
r1 = deg2km(r1);
filename = sprintf('xcorf/%s_%s_f',stainfo(sta1).staname,stainfo(sta2).staname);
data1 = load(filename);
xcorf1 = data1.cohere_sum./data1.coherenum;

r2 = distance(stainfo(sta3).lat,stainfo(sta3).lon,stainfo(sta4).lat,stainfo(sta4).lon);
r2 = deg2km(r2);
filename = sprintf('xcorf/%s_%s_f',stainfo(sta3).staname,stainfo(sta4).staname);
data2 = load(filename);
xcorf2 = data2.cohere_sum./data2.coherenum;

N = 1000;
xcorf1 = real(xcorf1(1:N));
xcorf2 = real(xcorf2(1:N));
%xcorf1 = smooth(xcorf1,10);
%xcorf2 = smooth(xcorf2,10);
faxis = [0:N-1]*1/3600;

waxis = (frange(1):1/3600:frange(2))*2*pi;

xsp1 = interp1(faxis*2*pi,xcorf1,waxis);
xsp2 = interp1(faxis*2*pi,xcorf2,waxis);

xsp1 = smooth(xsp1,10);
xsp2 = smooth(xsp2,10);

tw1 = ones(1,tN)*r1./refc;
tw1_0 = lsqnonlin(@(x) besselerr(x,[xsp1]),[tw1],[tw1]*0.8,[tw1]*1.2);
tw2 = ones(1,tN)*r2./refc;
tw2_0 = lsqnonlin(@(x) besselerr(x,[xsp2]),[tw2],[tw2]*0.8,[tw2]*1.2);

x1 = twloc.*tw1;
x2 = twloc.*tw2;
xmin = max([min(x1),min(x2)]);
xmax = min([max(x1),max(x2)]);
dx = (waxis(2)-waxis(1))*mean([tw1 tw2]);
xcommon = xmin:dx:xmax;

tw = lsqnonlin(@(x) xcorferr(x,[xsp1,xsp2]),[tw1_0 tw2_0],[tw1_0 tw2_0]*0.8,[tw1_0 tw2_0]*1.2);
tw1 = tw(1:tN);
tw2 = tw(tN+1:2*tN);

figure(2)
clf
hold on
plot(1./twloc*2*pi,r1./tw1,'ro-'); 
plot(1./twloc*2*pi,r2./tw2,'rx-'); 
plot(1./twloc*2*pi,r1./tw1_0,'ko-'); 
plot(1./twloc*2*pi,r2./tw2_0,'kx-'); 
% tw = lsqnonlin(@(x)xcorferr(x,[xsp1,xsp2]),[tw1 tw2]);


