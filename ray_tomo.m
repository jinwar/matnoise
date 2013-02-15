% Script to do the ray theory tomography based on the ambient noise measurement
% Written by Ge Jin, jinwar@gmail.com
% Nov 2012


% input files
load stainfo_BHZ.mat
load xspinfo.mat
load seiscmap.mat

% Set up geometry parameters
lalim=[-11.2 -7.8];
lolim=[148.8 151.5];
gridsize=0.1;
distrange= [2 6]; % in term of wavelength
smweight0 = 2;
maxerrweight = 2;
fiterrtol = 3;
errlevel = 1;
dterrtol = 2;
snrtol = 1.1;
isoutput = 1;
r=0.1;
refv = 3.2;
raydensetol=deg2km(gridsize)*2;

xnode=lalim(1):gridsize:lalim(2);
ynode=lolim(1):gridsize:lolim(2);
Nx = length(xnode);
Ny = length(ynode);
periods=2*pi./twloc;

% read in bad station list, if existed
if exist('badsta.lst')
	badstnms = textread('badsta.lst','%s');
	badstaids = find(ismember({stainfo.staname},badstnms));
	disp('Found Bad stations:')
	disp(badstnms)
end

% Set up initial smoothing kernel
[i,j] = ndgrid(1:Nx,2:(Ny-1));
ind = j(:) + Ny*(i(:)-1);
dy = diff(ynode);
dy1 = dy(j(:)-1);
dy2 = dy(j(:));
Areg = sparse(repmat(ind,1,3),[ind-1,ind,ind+1], ...
    [-2./(dy1.*(dy1+dy2)), 2./(dy1.*dy2), -2./(dy2.*(dy1+dy2))],Nx*Ny,Nx*Ny);
[i,j] = ndgrid(2:(Nx-1),1:Ny);
ind = j(:) + Ny*(i(:)-1);
dx = diff(xnode);
dx1 = dx(i(:)-1);
dx2 = dx(i(:));
Areg = [Areg;sparse(repmat(ind,1,3),[ind-Ny,ind,ind+Ny], ...
    [-2./(dx1.*(dx1+dx2)), 2./(dx1.*dx2), -2./(dx2.*(dx1+dx2))],Nx*Ny,Nx*Ny)];
F=Areg;

% Initial the xsp structure
for ixsp = 1:length(xspinfo)
	xspinfo(ixsp).isgood = 0;
	if xspinfo(ixsp).sumerr < errlevel ...
            && xspinfo(ixsp).snr > snrtol && xspinfo(ixsp).coherenum > 2000
        xspinfo(ixsp).isgood = 1;
	end
	if sum(ismember([xspinfo(ixsp).sta1 xspinfo(ixsp).sta2],badstaids)) > 0
		xspinfo(ixsp).isgood = 0;
	end
end % end of loop ixsp

for ip=1:length(periods)
	disp(' ');
	disp(['Inversing Period: ',num2str(periods(ip))]);
	clear rays dt fiterr mat phaseg err raydense
	raynum = 0;
	for ixsp = 1:length(xspinfo)
		if xspinfo(ixsp).isgood ==0
			continue;
		end
		if xspinfo(ixsp).r > refv*periods(ip)*distrange(2)...
				|| xspinfo(ixsp).r < refv*periods(ip)*distrange(1)
			continue;
		end
		raynum = raynum+1;
		rays(raynum,1) = stainfo(xspinfo(ixsp).sta1).lat;
		rays(raynum,2) = stainfo(xspinfo(ixsp).sta1).lon;
		rays(raynum,3) = stainfo(xspinfo(ixsp).sta2).lat;
		rays(raynum,4) = stainfo(xspinfo(ixsp).sta2).lon;
		dt(raynum) = xspinfo(ixsp).tw(ip);
		err = smooth((abs(xspinfo(ixsp).err)./mean(abs(xspinfo(ixsp).xsp))).^2,round(length(waxis)/length(twloc)));
		fiterr(raynum) = interp1(waxis(:),err(:),twloc(ip));
		csnum(raynum) = xspinfo(ixsp).coherenum;
        snr(raynum) = xspinfo(ixsp).snr;
	end
	if size(dt,1) ~=raynum
		dt = dt';
	end

	% Building the data kernel
	disp('Start building the kernel');
	tic
	mat=ray_kernel_build(rays,xnode,ynode);
	toc
	% Calculate the weighting matrix
	W = sparse(length(dt),length(dt));
	for i=1:length(dt)
		W(i,i)=1./fiterr(i);
	end
	ind = find(W > maxerrweight);
	W(ind) = maxerrweight;
	ind = find(W < 1/fiterrtol);
	W(ind) = 0;
	for i=1:length(dt)
		W(i,i)=W(i,i).*(csnum(i).^0.5);
	end

	% calculate the smoothing weight
	smweight = smweight0;
	NR=norm(F,1);
	NA=norm(W*mat,1);
	smweight = smweight0*NA/NR;
	
	disp('start inverse');
	A=[W*mat;smweight*F];
	rhs=[W*dt;zeros(size(F,1),1)];
	
%        disp('start inverse');
%        tic
	phaseg=(A'*A)\(A'*rhs);
%        toc
%        disp('Done');
	
	
	% Iteratively down weight the measurement with high error
	niter=1;
	
	while niter < 2
		niter=niter+1;
		err = mat*phaseg - dt;
%            err = W*err;
		stderr=std(err);
        if stderr > dterrtol
            stderr = dterrtol
        end
		ind = find(diag(W)==0);
		disp('Before iter:');
		disp(['Good Measurement Number: ', num2str(length(diag(W))-length(ind))]);
		disp(['Bad Measurement Number: ', num2str(length(ind))]);
		for i=1:length(err)
			if abs(err(i)) > 2*stderr
				W(i,i)=0;
			end
		end
		ind = find(diag(W)==0);
		disp('After iter:');
		disp(['Good Measurement Number: ', num2str(length(diag(W))-length(ind))]);
		disp(['Bad Measurement Number: ', num2str(length(ind))]);
		
		% Rescale the smooth kernel
		NR=norm(F,1);
		NA=norm(W*mat,1);
		smweight = smweight0*NA/NR;
		
		A=[W*mat;smweight*F];
		rhs=[W*dt;zeros(size(F,1),1)];
		
%            disp('start inverse');
%            tic
		phaseg=(A'*A)\(A'*rhs);
%            toc
%            disp('Done');
		
	end
	
%        disp(' Get rid of uncertainty area');
	for i=1:Nx
		for j=1:Ny
			n=Ny*(i-1)+j;
			raydense(i,j) = sum(mat(:,n));
			if raydense(i,j) < raydensetol
				phaseg(n)=NaN;
			end
		end
	end

	% Convert into phase velocity
	for i=1:Nx
		for j=1:Ny
			n=Ny*(i-1)+j;
			GV(i,j)= 1./phaseg(n);
		end
	end
	raytomo(ip).GV = GV;
	raytomo(ip).mat = mat;
	raytomo(ip).raydense = raydense;
	raytomo(ip).period = periods(ip);
end % end of period loop
lalim = [min(xnode) max(xnode)];
lolim = [min(ynode) max(ynode)];
[xi yi] = ndgrid(xnode,ynode);
if isoutput
save('raytomo.mat','raytomo','xnode','ynode');
save('coor.mat','xi','yi','xnode','ynode','gridsize','lalim','lolim');
end
figure(17)
clf
for ip=1:20
    subplot(4,5,ip)
    ax = worldmap(lalim, lolim);
    set(ax, 'Visible', 'off')
    surfacem(xi,yi,raytomo(ip).GV);
    drawpng
    title(['Periods: ',num2str(periods(ip))],'fontsize',15)
	avgv = nanmean(raytomo(ip).GV(:));
    caxis([avgv*(1-r) avgv*(1+r)])
    colorbar
	colormap(seiscmap)
end

figure(18)
clf
lalim = [min(xnode) max(xnode)];
lolim = [min(ynode) max(ynode)];
[xi yi] = ndgrid(xnode,ynode);
for ip=1:20
    subplot(4,5,ip)
    ax = worldmap(lalim, lolim);
    set(ax, 'Visible', 'off')
    surfacem(xi,yi,raytomo(ip).raydense);
    drawpng
    title(['Periods: ',num2str(periods(ip))],'fontsize',15)
    colorbar
end
