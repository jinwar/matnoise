function err = xcorferr(tw,xsps)
% Error function for fitting two nearby station pair's frequency domain xspectrum

global tN
global waxis
global twloc
global bzeros
global xcommon

sumerr = 0;
Isfigure=1;
interpmethod = 'linear';
tw1 = interp1(twloc,tw(1:tN),waxis,interpmethod,'extrap');
tw2 = interp1(twloc,tw(tN+1:2*tN),waxis,interpmethod,'extrap');

x1 = waxis.*tw1;
x2 = waxis.*tw2;

xsp1 = xsps(1:length(waxis));
xsp2 = xsps(length(waxis)+1:2*length(waxis));

% First part of error: two xspectrum should be like each other
F1 = interp1(x1,xsp1,xcommon,interpmethod,'extrap');
% F1 = F1./max(abs(F1));
F2 = interp1(x2,xsp2,xcommon,interpmethod,'extrap');
% F2 = F2./max(abs(F2));
be = besselj(0,xcommon);
be = be./mean(abs(be)).*mean([abs(F1) abs(F2)]);

% Second part of error: Fit Bessel Function

F1z =  be - F1; 
F2z =  be - F2; 

err = [F1-F2,F1z,F2z];

if Isfigure
    figure(1)
    clf
    hold on
    plot(xcommon,F1);
    plot(xcommon,F2,'r');
    plot(xcommon,be,'k');
end
disp(sum(err.^2))

end
