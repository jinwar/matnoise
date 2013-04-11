% setup parameters for ambient noise project

lalim=[-11.2 -7.8];
lolim=[148.8 151.5];
gridsize=0.1;
distrange= [2 6]; % in term of wavelength
smweight0 = 3;
maxerrweight = 2;
fiterrtol = 1;
errlevel = 1;
dterrtol = 3;
snrtol = 1.1;
polyfit_dt_err = 2;
isoutput = 1;
r=0.1;
refv = 3.2;
raydensetol=deg2km(gridsize)*2;

