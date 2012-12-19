sta1 = 1;
sta2 = 13;
sta3 = 1;
sta4 = 32;

r1 = distance(stainfo(sta1).lat,stainfo(sta1).lon,stainfo(sta2).lat,stainfo(sta2).lon);
filename = sprintf('xcorf/%s_%s_f',stainfo(sta1).staname,stainfo(sta2).staname);
data1 = load(filename);
xcorf1 = data1.cohere_sum./data1.coherenum;

r2 = distance(stainfo(sta3).lat,stainfo(sta3).lon,stainfo(sta4).lat,stainfo(sta4).lon);
filename = sprintf('xcorf/%s_%s_f',stainfo(sta3).staname,stainfo(sta4).staname);
data2 = load(filename);
xcorf2 = data2.cohere_sum./data2.coherenum;

N = 1800;
xcorf1 = xcorf1(1:N);
xcorf2 = xcorf2(1:N);
xcorf1 = smooth(xcorf1,10);
xcorf2 = smooth(xcorf2,10);
faxis = [0:N-1]*1/3600;

figure(1)
clf
subplot(2,1,1)
hold on
plot(faxis,abs(xcorf1));
plot(faxis,abs(xcorf2),'r');
subplot(2,1,2)
hold on
plot(faxis,real(xcorf1));
plot(faxis,real(xcorf2),'r');
