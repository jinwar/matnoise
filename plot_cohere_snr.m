clear;

load stainfo_BHZ.mat
load xspinfo.mat

snr = [xspinfo(:).snr];
r = [xspinfo(:).r];
coherenum = [xspinfo(:).coherenum];
sumerr = [xspinfo(:).sumerr];

figure(43);
clf
hist([xspinfo(:).snr],30);

figure(42)
clf
hist(coherenum)

figure(2);
clf
hold on
ind = find(r<100);
plot(coherenum(ind),snr(ind),'rx');
ind = find(r>100);
plot(coherenum(ind),snr(ind),'bx');
title('coherenum vs snr')

figure(3)
clf
hold on
plot(snr,sumerr,'x')
title('snr vs sumerr')

figure(4)
clf
ind = find(snr < 1.1 & sumerr<1);
length(ind)
ind = ind(end);

subplot(2,1,1)
plot(xspinfo(ind).xsp)
title(xspinfo(ind).filename)
subplot(2,1,2)
data1 =load(xspinfo(ind).filename);
xcorf1 = data1.cohere_sum./data1.coherenum;
faxis = [0:1799]/3600;
xcorf1 = xcorf1(1:length(faxis));
plot(faxis,real(xcorf1));
besselerr(xspinfo(ind).tw,xspinfo(ind).xsp,5);


