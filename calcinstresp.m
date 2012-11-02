function [resp dt]=calcinstresp(dbsnin,sta,chan,time,timelength, lo_corner)
%function [resp dt]=calcinstresp(dbsnin,sta,chan,time, fflen, lo_corner)
% From a station and channel and pointer to sensor-instrument join, get a
% response curve vs frequency
% input:  dbsnin = css3.0 db join of sensor - instrument
% (sta, chan) = station, channel
% time = time for which to sample the response (e.g. through inst changes)
%          <0 to just take the first instance of sta/chan
% fflen = # of points in output spectra; df = 1/(fflen*dt) 
%       NOTE  gets dt from sensor table
% lo_corner = Low frequency (high-pass) corner frequency for output
% sets response to 1 in velocity at >0.5*(nyquist_frequency) to leave FIR
% antialias filters alone
%   OUTPUT resp is instrument response per unit ground DISPLACEMENT
%     expressed as transfer function:  (timeseries)*(resp)=(ground displacement)
% GA 3/09
%
% Modefied by Ge Jin, Oct 2012 to change the input from points of data to time length
bufnyq=0.5;
npoles=5;    % Nominal filter cutoff for response in Butterworth poles (acausal)
if (time<0)
   dbsn=dbsubset(dbsnin,sprintf('sta=~/%s/ && chan=~/%s/',sta,chan));
else
   dbsn=dbsubset(dbsnin,sprintf('sta=~/%s/ && chan=~/%s/ && time<%.3f && endtime>%.3f ',...
       sta,chan,time,time));
end
dbsn.record=0;   % Kluge. fix later:  should hopefully just have one
instfile=dbfilename(dbsn);
respobj=dbresponse(instfile);
    % or dbextfile on join w/ wfdisc?
[ncalib rsptype samprate] = dbgetv(dbsn,'ncalib','rsptype','samprate');
dt=1./samprate;
fflen = round(timelength*samprate);
jny1=fflen/2;
frq=2.*pi.*(0:jny1)./(fflen.*dt);
frq=[frq,-fliplr(frq(2:jny1))]';     % adjust to FFT convention
knyq=find(abs(frq)>bufnyq*frq(jny1+1));   % flag frequencies more than twice the nyquist

if ncalib<0, ncalib=1; end;
if (rsptype=='D' || rsptype=='d')
    acor=ones(size(frq))./ncalib;
elseif (rsptype=='A')
    acor=1./(i.*frq)./ncalib;
else
    acor=i.*frq./ncalib;    % corrects nominal velocity response to displacement
end
resp=eval_response(respobj,frq).*acor;   % complex response fcns over fftlen...
resp(knyq)=acor(knyq);
lo_w=2*pi*lo_corner;
hpfiltfrq=( ((frq./lo_w).^(2*npoles))./(1+(frq./lo_w).^(2*npoles)) );
resp=hpfiltfrq./resp;    % this is normalization transfer function
resp(1)=0;
knan=find(isnan(resp));
if ~isempty(knan), resp(knan)=0; end;
return
