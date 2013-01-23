% This script is used to calculate the BHD to BHZ transfer function to correct the compliance.
% Method is descripted in the paper: Crawford and Webbs, 2000
% written by Ge Jin

clear;

TimeSegN = 1000;
load stainfo_BHZ.mat

for ista = 36:43   % index of OBSs
	staname = stainfo(ista).staname;
	disp(staname)
	Nd = 0;
	for itime = 0:11
		filename = ['data/',staname,'_',num2str(itime),'.mat'];
		bhz = load(filename);
		filename = ['data/',staname,'_',num2str(itime),'_BDH.mat'];
		bdh = load(filename);
		for iseg = itime*TimeSegN:(itime+1)*TimeSegN-1
			if iseg == 0
				continue;
			end
			if stainfo(ista).datacover(iseg) == 1
				Nd = Nd+1;
				data_z = bhz.segdata(iseg,:);
				data_h = bdh.segdata(iseg,:);

				data_z_f = fft(data_z);
				data_h = fft(data_z);


			end
		
		end % end of iseg
	end % end of itime
end  %end of stations

