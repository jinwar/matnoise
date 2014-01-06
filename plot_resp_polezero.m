% matlab script to make the plot of instrument response based on poles and zeros

clear;
freq = [0:0.0001:0.1];

%station 1
gain = 5.730e8;
i = sqrt(-1);
poles = [-3.701e-2 + i*3.701e-2;...
		 -3.701e-2 - i*3.701e-2;...
	-1.131E+03; -1.005E+03;-5.027E+02 ];
zeros = [0;0];

w = freq.*2*pi;
resp = ones(size(w));
		
for ip = 1:length(poles)
	resp = resp./(i*w - poles(ip));
end
for ip = 1:length(zeros)
	resp = resp.*(i*w - zeros(ip));
end
resp = resp*gain;

figure(25)
clf
hold on
subplot(1,2,1)
hold on
plot(freq,abs(resp),'rx');
subplot(1,2,2)
hold on
plot(freq,angle(resp),'rx');

% station 2
gain = 4.532E+05;
i = sqrt(-1);
poles = [ ...
-0.01813   +i*  0.01803;...
-0.01813   +i* -0.01803;...
-124.9     +i*  0.0000 ;...
-197.5     +i*  256.1  ;...
-197.5     +i* -256.1  ;...
-569       +i*  1150   ;...
-569       +i* -1150   ;...
]
zeros = [0;0;-90;-164.2;-3203];
w = freq.*2*pi;
resp = ones(size(w));
		
for ip = 1:length(poles)
	resp = resp./(i*w - poles(ip));
end
for ip = 1:length(zeros)
	resp = resp.*(i*w - zeros(ip));
end
resp = resp*gain;

figure(25)
subplot(1,2,1)
hold on
plot(freq,abs(resp),'bx');
subplot(1,2,2)
hold on
plot(freq,angle(resp),'bx');


