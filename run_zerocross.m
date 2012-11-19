% Script to test the bessel zero crossing method
% function is coded by Yang Zha
clear
load ref_ocean1500.mat
load xspinfo.mat
load stainfo_BHZ.mat
load refphasev.mat

warning off

ref_dispersion(:,1) = f_ref  ;
ref_dispersion(:,2) = c_ref  ;
ref_dispersion(:,3) = vg_ref ;
ref_dispersion(:,4) = grad_c ;
ref_dispersion(:,5) = grad_vg;


figure(1)
clf
hold on

for ixp = 1:length(xspinfo)
	disp(ixp)
	if xspinfo(ixp).r > 80
        sta1 = xspinfo(ixp).sta1;
        sta2 = xspinfo(ixp).sta2;
        filename = sprintf('xcorf/%s_%s_f.mat',stainfo(sta1).staname,stainfo(sta2).staname);
        if ~exist(filename,'file')
            disp(['not exist ',filename])
            continue;
        end
        data1 = load(filename);
        xcorf1 = data1.cohere_sum./data1.coherenum;
		coh(:,1) = [0:1799]*1/3600;
		coh(:,2) = smooth(real(xcorf1(1:1800)),10);
		[freq, c, vg] = CalcPhaseV_Aki([0.04 0.25],xspinfo(ixp).r,ref_dispersion,coh);
%		plot(freq,c,'-o');
%        drawnow;
        if sum(isnan(freq(:)))==0;
            xspinfo(ixp).zctw = interp1(freq,xspinfo(ixp).r./c,twloc/2/pi);
        end
    end
end
