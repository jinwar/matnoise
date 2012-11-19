function [freq, c, vg] = CalcPhaseV_Aki(freq_range, stadist, ref_dispersion, coh, varargin)

%% A function to recover phase and group velocity from coherence of ambient
% noise data using AKI's formulism 

% INPUT: 
% freq_range: range of frequency to pick, [f_min, f_max].
% stadist: station distance in km
% ref_dispersion: a reference dispersion curve(both vg and vc with their
% respective slope), to set a boundary for phase velocity picking.
% ref_dispersion: [f_ref,c_ref,vg_ref,grad_c, grad_vg]
% coh : smoothed coherence function and frequency.
% options: 

% OUTPUT:
% c: phase velocity in km/s
% vg: group velocity
% Freq: frequency at which c and vg are given

%% setup default parameters
vec_m=[-2:1:2];
IFplot = false;
vg_range=[1,1]; 
% [lower upper] range of allowed group velocity variation from reference model
sigma = 10;% slope_parameter, tolerance=range*sqrt(1+slope^2*sigma) (km/s); 
% when sigma=0; tolerance of velocity range is same for all frequencies. 
c_range=[0.7,0.7];
f0= 1/5;
%% input argument

if nargin==5 % if option is provided
     
     option = varargin{1}; 
     vg_range = option.vg_range;
     c_range = option.c_range;
     vec_m = option.vec_m;
     sigma = option.sigma;
     IFplot = option.IFplot;

end

 
%% Pick zero crossings in coherence-frequency relation


vec_f = coh(:,1);
coherence = coh(:,2); % smoothed coherence;


%range=[1.5*v_min/stadist 0.5];

%ave_coh_rl = fastsmooth(coh(:,2),navg,2); % moving average every 100
%points



% find zeros and calcualte phase velocity at corresponding freqency
[ind,freqzero,nzero,sign]=FindArrayZero(vec_f, coherence, freq_range);

% if the first crossing is from negative to positive, then discard it.
if sign(1)==1 
    freqzero=freqzero(2:length(freqzero));
        sign=sign(2:length(sign));
        nzero=nzero-1;
end

t_zero=1./freqzero; % period 

%% Formulate a suite of dispersion curves

mtx_c = zeros(nzero,length(vec_m));

for im=1:length(vec_m)
c_phase = zeros(nzero,1);
m=vec_m(im);

if m<0
    jzeros=besselzero(0,nzero+2*m,1);
    shift_jzeros=[ones(-2*m,1)*nan;jzeros(1:nzero+2*m)];
    c_phase=stadist*2*pi*freqzero./shift_jzeros;
end

if m>=0
    
    jzeros=besselzero(0,nzero+2*m,1);
    shift_jzeros=jzeros(2*m+1:nzero+2*m);

    c_phase=stadist*2*pi*freqzero./shift_jzeros; % phase velocity in km/s
end

% average up and down crossings
c_up=interp1(freqzero(sign==1),c_phase(sign==1),freqzero);
c_dn=interp1(freqzero(sign==-1),c_phase(sign==-1),freqzero);
mtx_c(:,im)=0.5 * (c_up+c_dn); % average up abd down crossings
lgd{im}=strcat('m= ',num2str(m));
end



%% Pick the phase velocity curve based on the following criteria
% 1. the phase velocity increase with period 
% 2. the phase velocity and the calculated group velocity fall into the
% range provided by the reference model
% 3. if 1 and 2 are satisfied, choose the curve with minimum local
% curvature (minimum difference in slope);

% ------------------------------------------------------------------------
% load reference curve
f_ref = ref_dispersion(:,1);
c_ref = ref_dispersion(:,2);
vg_ref = ref_dispersion(:,3);
grad_c = ref_dispersion(:,4);
grad_vg = ref_dispersion(:,5);

v_ref=interp1(f_ref,c_ref,freqzero,'linear','extrap');
c_slope=interp1(f_ref,grad_c,freqzero,'linear','extrap');



n=1;% index of valid frequency to evaluate dispersion curve
for i=1:length(freqzero)
    % set upper and lower bound of allowed velocity
    vmin(i)=max(v_ref(i)-c_range(1)*sqrt(1+sigma*c_slope(i)^2), 1.5);
    vmax(i)=v_ref(i)+c_range(2)*sqrt(1+sigma*c_slope(i)^2);
    

    switch n
        case 1 % if this is the first valid velocity to be determined
           condition=((mtx_c(i,:)>vmin(i)).*(mtx_c(i,:)<vmax(i))); % only apply velocity bound filter 
           if sum(condition)==1
               c(n)=mtx_c(i,find(condition));
               freq(n)=freqzero(i);
               k(n)=2*pi*freq(n)/c(n);
               n=n+1;
               
           end
        case 2       
           condition1 = ( mtx_c(i,:)<c(n-1)*1.05 ); % velocity decrease with increasing frequency 
           
           condition2=((mtx_c(i,:)>vmin(i)).*(mtx_c(i,:)<vmax(i))); % velocity bound 
           k_number=2*pi*freqzero(i)./mtx_c(i,:); % trial wave number
               
               vec_vg=2*pi*(freqzero(i)-freq(n-1))./(k_number-k(n-1)); % measured group velocity
               
               vg_slope=interp1(f_ref,grad_vg,0.5*(freqzero(i)+freq(n-1)),'linear','extrap');
               
               vg_predict=interp1(f_ref,vg_ref,0.5*(freqzero(i)+freq(n-1)),'linear','extrap'); 
               
               vg_max(i)=vg_predict+vg_range(2)*sqrt(1+sigma*vg_slope^2);
               vg_min(i)=vg_predict-vg_range(1)*sqrt(1+sigma*vg_slope^2);
               
               condition3=(vec_vg<vg_max(i)).*(vec_vg>vg_min(i));
               
           condition = condition1.*condition2.*condition3;
           %    phase velocity must decrease with frequency, with a little
           %    tolerance v(n)<v(n-1)*1.05
           if sum(condition)>=1
               vec_c=mtx_c(i,find(condition));
               condition3 = abs(vec_c-c(n-1))==min(abs(vec_c-c(n-1))); % minimal decrease each step. 
               sum(condition3);
               c(n)=vec_c(find(condition3));
               freq(n)=freqzero(i);
               k(n)=2*pi*freq(n)/c(n);
               n=n+1;
           end
        otherwise
               condition1 = ( mtx_c(i,:)<c(n-1)*1.05 );
               condition1b = ( mtx_c(i,:)<c(n-1) );
               condition2=((mtx_c(i,:)>vmin(i)).*(mtx_c(i,:)<vmax(i))); % phase velocity within range
               condition4 = (mtx_c(i,:)>1.5) ;% phase velocity > that of water
               % condition3: predicted group velocuty must be within range
               
               
               k_number=2*pi*freqzero(i)./mtx_c(i,:); % trail wave number
               
               vec_vg=2*pi*(freqzero(i)-freq(n-1))./(k_number-k(n-1)); % measured group velocity
               
               vg_slope=interp1(f_ref,grad_vg,0.5*(freqzero(i)+freq(n-1)),'linear','extrap');
               
               vg_predict=interp1(f_ref,vg_ref,0.5*(freqzero(i)+freq(n-1)),'linear','extrap'); 
               
               vg_max(i)=vg_predict+vg_range(2)*sqrt(1+sigma*vg_slope^2);
               vg_min(i)=vg_predict-vg_range(1)*sqrt(1+sigma*vg_slope^2);
               
               condition3=(vec_vg<vg_max(i)).*(vec_vg>vg_min(i));% group velocity in range
               if freqzero(i) < f0
               condition = condition1.*condition3.*condition4; 
               else
                   condition = condition1.*condition2.*condition3.*condition4;
               end
               % predicted vg must be within +/- 1km/s of reference model
               
               if sum(condition)>=1 % if all the above conditions are met, choose a point that yield a smoothest curve
                   vec_c=mtx_c(i,find(condition));
                   slope1 = ( c(n-1) - c(n-2) ) / ( freq(n-1) - freq(n-2) );
                   slope2 = (vec_c - c(n-1)) / (freqzero(i)-freq(n-1));
                   condition3 = abs(slope2-slope1)==min(abs(slope2-slope1)); % minimal decrease each step. 
                   c(n)=vec_c(find(condition3));
                   freq(n)=freqzero(i);
                   k(n)=2*pi*freq(n)/c(n);
                   n=n+1;
                end
     end
           

   
end

if n>3
    

%% Plot phase and group velocity
n
vec_k=2*pi*freq./c;
freq_vg = 0.5* ( freq(1:length(freq)-1)+freq( 2:length(freq)) );
vg= 2*pi*diff(freq)./diff(vec_k);
vg = interp1(freq_vg,vg,freq,'linear','extrap');

if(IFplot)
    
figure;hold on

plot(t_zero,mtx_c,'-o','linewidth',1);


%legend('m=0','m=1','m=2','m=-1','m=-2')
legend(lgd)
%xlim([1,25]);
ylim([0.5,5]);
%title(strcat(sta1,'-',sta2));
%plot(freqzero,vmax,'-.r','linewidth',0.5);
%plot(freqzero,vmin,'-.r','linewidth',0.5);    
plot(1./freq, c ,'-r*','linewidth',2,'markersize',10);
plot(1./freqzero,vmax,'-.r','linewidth',0.5);
plot(1./freqzero,vmin,'-.r','linewidth',0.5);
%plot(f_ref,vg_ref+vg_range(2).*sqrt(1+sigma*grad_vg.^2),'-.b','linewidth',0.5);
%plot(f_ref,vg_ref-vg_range(1).*sqrt(1+sigma*grad_vg.^2),'-.b','linewidth',0.5);
%xlim([1,40]);
xlabel('frequency / Hz')
ylabel('Phase velocity / km/s')


%figure; hold on;
%plot(1./freq,c ,'-r*','linewidth',2,'markersize',5);
%plot(1./freq,vg,'-b*','linewidth',2,'markersize',5);
%legend('Phase velocity','group velocity');
%% compare with reference model
%
%plot(1./freqzero,vmax,'-.r','linewidth',0.5);
%plot(1./freqzero,vmin,'-.r','linewidth',0.5);
%plot(1./f_ref,vg_ref+vg_range(2).*sqrt(1+sigma*grad_vg.^2),'-.b','linewidth',0.5);
%plot(1./f_ref,vg_ref-vg_range(1).*sqrt(1+sigma*grad_vg.^2),'-.b','linewidth',0.5);
%ylim([0.5,5]);xlim([2,25]);
%fprintf('station distance is %f km \n',stadist);
%%title(strcat(sta1,'-',sta2));
%
%
%    
%figure; hold on;
%subplot(3,2,1)
%plot(freq,c,'-r*','linewidth',1,'markersize',5);
%title('phase_velocity');
%xlabel('frequency / Hz');
%
%subplot(3,2,2)
%plot(1./freq,c,'-r*','linewidth',1,'markersize',5);
%title('phase_velocity');
%xlabel('Period /s ');
%
%subplot(3,2,3);
%
%plot(freq,vec_k,'-ro','linewidth',1,'markersize',5);
%title('wave number k');
%xlabel('frequency / Hz');
%
%subplot(3,2,4)
%plot(1./freq,vec_k,'-r*','linewidth',1,'markersize',5);
%title('wave numiber');
%xlabel('Period /s ');
%
%
%
%subplot(3,2,5);
%
%plot(freq,vg,'-ro','linewidth',1,'markersize',5);
%title('group velocity dw/dk');
%xlabel('frequency / Hz');
%
%subplot(3,2,6)
%plot(1./freq,vg,'-r*','linewidth',1,'markersize',5);
%title('group velocity');
%xlabel('Period /s ');


end
else
freq=nan*ones(size(f_ref));
    c=freq;
    vg=c;
end


freq=freq(:);
c=c(:);
vg=vg(:);
return 
end
