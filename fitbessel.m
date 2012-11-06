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
refc2 = 2.5;
refc = 3.5 + (refc2-refc1)/(tN)*[1:tN];

twloc = frange(1):(frange(2)-frange(1))/(tN-1):frange(2);
twloc = twloc*2*pi;
bzeros = besselzero(0,100,1);

sta1 = 1;
sta2 = 11;
r1 = distance(stainfo(sta1).lat,stainfo(sta1).lon,stainfo(sta2).lat,stainfo(sta2).lon);
r1 = deg2km(r1);
filename = sprintf('xcorf/%s_%s_f',stainfo(sta1).staname,stainfo(sta2).staname);
data1 = load(filename);
xcorf1 = data1.cohere_sum./data1.coherenum;

N = 1000;
xcorf1 = real(xcorf1(1:N));
%xcorf1 = smooth(xcorf1,10);
%xcorf2 = smooth(xcorf2,10);
faxis = [0:N-1]*1/3600;

waxis = (frange(1):1/3600:frange(2))*2*pi;

xsp1 = interp1(faxis*2*pi,xcorf1,waxis);

xsp1 = smooth(xsp1,10);

tw1 = ones(1,tN)*r1./refc;

x1 = twloc.*tw1;
xmin = min(x1);
xmax = max(x1);
dx = (waxis(2)-waxis(1))*mean([tw1]);
xcommon = xmin:dx:xmax;

tw = lsqnonlin(@(x) besselerr(x,[xsp1]),[tw1],[tw1]*0.8,[tw1]*1.2);
% tw = lsqnonlin(@(x)xcorferr(x,[xsp1,xsp2]),[tw1 tw2]);
figure(2)
clf
plot(1./twloc*2*pi,r1./tw);

