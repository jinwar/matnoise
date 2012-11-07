function err = besselerr(tw,xsp,varargin)
% Error function for fitting two nearby station pair's frequency domain xspectrum

global tN
global waxis
global twloc
global weight 

Isfigure=0;
interpmethod = 'linear';


if nargin>2 % if option is provided
     Isfigure = varargin{1};
end

tw1 = interp1(twloc,tw(1:tN),waxis,interpmethod);

x1 = waxis.*tw1;

% First part of error: fit the Bessel Function
F1 = xsp;
be = besselj(0,x1);
be = be./mean(abs(be)).*mean([abs(F1)]);
F1z =  be(:) - F1(:); 
F1z = F1z.*weight(:);
% Second part of error: dispersion curve smoothness
sm = del2(tw1);
sm = sm./mean(abs(sm)).*mean(abs(F1z));

err = [F1z(:); sm(:)*0.2];

% err = err./mean(abs(err))*1;

if Isfigure>0
    figure(Isfigure)
    clf
    subplot(2,1,1)
    hold on
    plot(waxis/2/pi,F1);
    plot(waxis/2/pi,be,'r');
	title( num2str(sum(err.^2)./length(err)./sum(F1.^2)*length(F1)));
    subplot(2,1,2)
    plot(waxis/2/pi,F1z(:).^2,'k')
end
end