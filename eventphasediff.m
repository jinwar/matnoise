% This code is used to convert the result of ambient noise frequency method bessel function
% fitting phase measurement into a event based phase difference measurement between nearby stations
% as the input for the Eikonal Tomography 
% Written by Ge Jin
% jinwar@gmail.com

clear

load stainfo_BHZ.mat
load xspinfo.mat

% set good bessel fitting
errlevel = 1;
nbdist = 150;
refv =  3.0;
Isfigure = 1;

for ixsp = 1:length(xspinfo)
	xspinfo(ixsp).isgood = 0;
	if xspinfo(ixsp).sumerr < errlevel
        dx = diff(xspinfo(ixsp).tw.*twloc);
        badind = find(dx<0);
        if length(badind)==0
			xspinfo(ixsp).isgood = 1;
		end
	end
end % end of loop ixsp

for ista = 1:length(stainfo)
% for ista = 4
	disp(stainfo(ista).staname);
	evla = stainfo(ista).lat;
	evlo = stainfo(ista).lon;
	xspind = find([xspinfo(:).sta1] == ista | [xspinfo(:).sta2] == ista );
	evid = ista;


	% set up station array
	stanum = 0;
	for ixsp = xspind
		if xspinfo(ixsp).isgood == 0
			continue;
		end
		if xspinfo(ixsp).sta2 == ista
			staid = xspinfo(ixsp).sta1;
		else
			staid = xspinfo(ixsp).sta2;
		end
		stanum = stanum + 1;
		slat(stanum) = stainfo(staid).lat;
		slon(stanum) = stainfo(staid).lon;
		staids(stanum) = staid;
		xspids(stanum) = ixsp;
		epidist(stanum) = deg2km(distance(evla,evlo,slat(stanum),slon(stanum)));
	end % end of xsp loop

	% Loop through station array and find nearby connections
	csnum=0;
    clear ray dt ddist fiterr
	for ista1 = 1:stanum
		sta1_id = staids(ista1);
		tw1 = xspinfo(xspids(ista1)).tw;
		err1 = xspinfo(xspids(ista1)).err;
		xsp1 = xspinfo(xspids(ista1)).xsp;
		dist = deg2km(distance(slat,slon,slat(ista1),slon(ista1)));
		nearstaind = find(dist<nbdist & dist > 1);
		for ista2 = nearstaind
			sta2_id = staids(ista2);
			if sta1_id > sta2_id
				continue;
			end
			tw2 = xspinfo(xspids(ista2)).tw;
			err2 = xspinfo(xspids(ista2)).err;
			xsp2 = xspinfo(xspids(ista2)).xsp;
			csnum = csnum+1;
			ray(csnum,1) = slat(ista1);
			ray(csnum,2) = slon(ista1);
			ray(csnum,3) = slat(ista2);
			ray(csnum,4) = slon(ista2);
			dt(csnum,:) = tw1(:) - tw2(:);
			ddist(csnum) = epidist(ista1) - epidist(ista2);
			fiterr(csnum,:) = (err1(:)./mean(abs(xsp1))).^2 + (err2(:)./mean(abs(xsp2))).^2;


		end %end of ista2 loop
	end % end of ista1 loop

    % add in station pairs from main station if the epidist is smaller than nbdist
	nearstaind = find(epidist < nbdist);
	for ista1 = nearstaind
		tw1 = xspinfo(xspids(ista1)).tw;
		err1 = xspinfo(xspids(ista1)).err;
		xsp1 = xspinfo(xspids(ista1)).xsp;

		csnum = csnum+1;
		ray(csnum,1) = slat(ista1);
		ray(csnum,2) = slon(ista1);
		ray(csnum,3) = evla;
		ray(csnum,4) = evlo;
		dt(csnum,:) = tw1;
		ddist(csnum) = epidist(ista1);
		fiterr(csnum,:) = 2*(err1(:)./mean(abs(xsp1))).^2;
	end % end of main station nbdist
    
    % correct cycle skipping
    for ics = 1:csnum
        syndt = ddist(ics)./refv;
        N = -2:2;
        for iw=1:length(tw1)
            testdt = dt(ics,iw) + N.*2*pi./twloc(iw);
            dterr = abs(testdt - syndt);
            [bestdterr bestdti] = min(dterr);
            dt(ics,iw) = testdt(bestdti);
            bestcycle(ics,iw) = bestdti;
        end % end of periods loop
    end

    event(ista).ray = ray;
    event(ista).dt = dt;
    event(ista).ddist = ddist;
    event(ista).fiterr = fiterr;
    event(ista).bestcycle = bestcycle;
	
    if Isfigure && csnum > 0
        plotcsray
        pause(0.5)
    end

end % end of ista loop

save('events.mat','event');
