function err = besselerr(tw,xsp)
% Error function for fitting two nearby station pair's frequency domain xspectrum

global tN
global waxis
global twloc
global bzeros
global xcommon

sumerr = 0;
Isfigure=1;
interpmethod = 'spline';
tw1 = interp1(twloc,tw(1:tN),waxis,interpmethod);

x1 = waxis.*tw1;

xsp1 = xsp(1:length(waxis));

% First part of error: two xspectrum should be like each other
F1 = interp1(x1,xsp1,xcommon,interpmethod);
% F1 = F1./max(abs(F1));
be = besselj(0,xcommon);
be = be./mean(abs(be)).*mean([abs(F1)]);
% Second part of error: Fit Bessel Function
F1z =  be - F1; 
err = [F1z];

if Isfigure
    figure(1)
    clf
    hold on
    plot(xcommon,F1);
    plot(xcommon,be,'k');
	title( num2str(sum(err.^2)));
end


end
