% this code is used to compare the coherency between BDH and BHZ component to see whether it's reasonable to correct
% the long frequency noise in BHZ component by using BDH data
clear;

load stainfo_BHZ.mat

staid = 43
staname = stainfo(staid).staname

timeseg = 4;

filename = ['data/',staname,'_',num2str(timeseg),'.mat'];
bhz = load(filename);

filename = ['data/',staname,'_',num2str(timeseg),'_BDH.mat'];
bdh = load(filename);

nhour = 1:400;
data1 = bhz.segdata((timeseg)*1000:timeseg*1000+999,:);
data1 = data1(nhour,:)';
data1 = data1(:);
data1 = detrend(data1);
data1 = data1./mean(abs(data1));

data2 = bdh.segdata((timeseg)*1000:timeseg*1000+999,:);
data2 = data2(nhour,:)';
data2 = data2(:);
data2 = detrend(data2);
data2 = data2./mean(abs(data2));

fN = 0.5;
T1 = 10;
T2 = 100;
w2 = 1/T1/fN;
w1 = 1/T2/fN;
[b a] = butter(2,[w1 w2]);
fdata1 = filtfilt(b,a,data1);
fdata2 = filtfilt(b,a,data2);

faxis = [0:floor(length(data1)/2)]/length(data1);
spectrum1 = fft(data1);
spectrum1 = abs(spectrum1(1:length(faxis)));
spectrum2 = fft(data2);
spectrum2 = abs(spectrum2(1:length(faxis)));
%
%
disp('Calculating the coherency');
winlength = 3600;
[cxy, f] = mscohere(data1,data2,winlength);
f = f/2/pi;

figure(1)
clf
plot(f,cxy);
set(gca,'fontsize',20)
xlabel('Frequency (Hz)');
ylabel('Coherency ');
title(staname);

figure(2)
clf
hold on
plot(fdata1)
plot(fdata2,'r')
xlim([100 5000])

figure(3)
clf
hold on
plot(data1)
plot(data2,'r')
xlim([100 5000])
%
figure(4)
clf
semilogy(faxis,spectrum1);
hold on
semilogy(faxis,spectrum2,'r');

