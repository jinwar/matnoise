% use a 2 layer crustal + mantle model to generate a reference phase velocity curve for use by CalcPhaseV_Aki
%
clear
load xspinfo.mat
load refphasev.mat

vec_T=[2:25];
%h_mantle=20.0;
vpvs = 1.9;
h_crust=30.0;
%depth=[0;cumsum(refmod1(:,1))];
vec_h = [5 20 20];
vec_vp = [5 6.4 8.0];
vec_vs = vec_vp/vpvs;
vec_rho = [2.7 2.8 3.2];

refmod(:,1) = vec_h(:);
refmod(:,2) = vec_vp(:);
refmod(:,3) = vec_vs(:);
refmod(:,4) = vec_rho(:);

% write model to file
%writemod_surf96(refmod,'ref.mod');
%c_surf96 = dispR_surf96(vec_T,refmod);

[c, vg, grad_c, grad_vg] = Calc_Ray_dispersion(vec_T,refmod,1,0);

figure(54);
clf
hold on
plot(vec_T,c,'-ro');
plot(2*pi./twloc,refv,'-bo');

xlabel('period / seconds');
ylabel('phase velocity / km/s')
xlim([6 25])

save 3layermodel.mat refmod c vg grad_c grad_vg vec_T
