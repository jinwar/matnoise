% This script is to read in the event structure and calculate the eikonal tomography for each event
% Ge Jin, jinwar@gmail.com
clear

% Input event list file
load events.mat
load stainfo_BHZ.mat
load xspinfo.mat
load refphasev.mat

% some constants
ERRTOR=0.5;			 % the error allowed for cs measurement
mincsnum=100;
sou_dist_tol = 1;  % count by wavelength
isfigure=0;

%phvrange(1,:)=[3.55 4.15];
periods=2*pi./twloc;

smweight0 = 1;
maxerrweight =2;

lalim=[-11.2 -7.8];
lolim=[148.8 151.5];
gridsize=0.1;

%lalim=[30 50];
%lolim=[-125 -90];
%gridsize=0.3;

raydensetol=deg2km(gridsize)*1;
lolim=lolim;
lat0=mean(lalim);
lon0=mean(lolim);
xnode=lalim(1):gridsize:lalim(2);
ynode=lolim(1):gridsize:lolim(2);

Nx=length(xnode);
Ny=length(ynode);
[xi yi]=ndgrid(xnode,ynode);

% new, using Laplacian term
disp('initial the smoothing kernel')
tic
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

F=sparse(Nx*Ny*2*2,Nx*Ny*2);
oldprocess=0;
for n=1:size(Areg,1)
%    newprocess=floor(n/2/Nx/Ny*10);
%    if newprocess>oldprocess
%        disp(['Finished: ', num2str(newprocess), '0%'])
%        oldprocess=newprocess;
%    end
    ind=find(Areg(n,:)~=0);
    F(2*n-1,2*ind-1)=Areg(n,ind);
    F(2*n,2*ind)=Areg(n,ind);
end
toc



for ie = 1:size(event,1)
%for ie = 1
    for ip=1:length(periods)
        
        disp(['Event ID: ',num2str(ie)]);
        disp(['Period: ',num2str(periods(ip))]);
        
        %lalim=[27 50.5];
        
        
        % read in data and information
        
        csnum=event(ie,ip).csnum;
        if csnum < mincsnum
			event_tomo(ie,ip).GV = zeros(size(xi));
			event_tomo(ie,ip).GV(:) = NaN;
			event_tomo(ie,ip).raydense = zeros(size(xi));
            continue;
        end
        
        % Make the matrix
        %mat=sparse(csnum,stanum);
        dt=event(ie,ip).dt;
		if size(dt,1) == 1
			dt = dt';
		end
        rays=event(ie,ip).ray;
        W = sparse(length(dt),length(dt));
        for i=1:length(dt)
            W(i,i)=1./event(ie,ip).fiterr(i);
        end
        ind = find(W > maxerrweight);
        W(ind) = maxerrweight;
        
        for i=1:csnum
            
            Isinmap=1;
            temp=[rays(i,1) rays(i,3)];
            if min(temp) < lalim(1) || max(temp) > lalim(2)
                Isinmap=0;
            end
            temp=[rays(i,2) rays(i,4)];
            if min(temp) < lolim(1) || max(temp) > lolim(2)
                Isinmap=0;
            end
            if ~Isinmap
                rays(i,1:4)=0;
                dt(i)=0;
                W(i,i)=0;
            end
        end
        disp('Start building the kernel');
        tic
        mat=kernel_build(rays,xnode,ynode);
        toc
        
        % Calculate the kernel density
        %sumG=sum(abs(mat),1);
        ind=1:Nx*Ny;
        sumG(ind)=sum((mat(:,2*ind).^2+mat(:,2*ind-1).^2).^.5,1);
		clear raydense
        for i=1:Nx
            for j=1:Ny
                n=Ny*(i-1)+j;
                raydense(i,j)=sumG(n);
            end
        end
        
        % build the smoothing operator
        
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
            niter
            niter=niter+1;
            err = mat*phaseg - dt;
%            err = W*err;
            stderr=std(err);
            for i=1:length(err)
                if abs(err(i)) > 2*stderr
                    W(i,i)=0;
                end
            end
			ind = find(diag(W)==0);
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
                if raydense(i,j) < raydensetol
                    phaseg(2*n-1)=NaN;
                    phaseg(2*n)=NaN;
                end
            end
        end
        
        % Change phaseg into phase velocity
        for i=1:Nx
            for j=1:Ny
                n=Ny*(i-1)+j;
                GVx(i,j)= phaseg(2*n-1);
                GVy(i,j)= phaseg(2*n);
            end
        end
        
        GV=(GVx.^2+GVy.^2).^-.5;

%       Get rid of the area that is too close to the source
		dist = deg2km(distance(xi,yi,event(ie,ip).evla,event(ie,ip).evlo));
		ind = find(dist < sou_dist_tol*periods(ip)*refv(ip));
		GV(ind) = NaN;

		event_tomo(ie,ip).GV = full(GV);
		event_tomo(ie,ip).raydense = full(raydense);
        
    end	% end of period loop
end % end of event loop
save('event_tomo.mat','event_tomo','xnode','ynode','periods');
