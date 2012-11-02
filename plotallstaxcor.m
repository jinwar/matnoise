clear

load stainfo_BHZ
for ista = 1:length(stainfo)
    plotsinglestaxcor(ista);
    figure(1)
    print('-dpng',stainfo(ista).staname);
end