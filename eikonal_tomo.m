% This script is to read in the event structure and calculate the eikonal tomography for each event
% Ge Jin, jinwar@gmail.com
clear

% Input event list file
load events.mat
load stainfo_BHZ.mat
load xspinfo.mat
load refphasev.mat
load raytomo.mat

% some constants
mincsnum=50;
sou_dist_tol = 1;  % count by wavelength
smweight0 =1.0;
Rdumpweight0 = 0.1;
Tdumpweight0 = 0.2;
maxerrweight = 2;
fiterrtol = 3;
dterrtol = 1;
isRsmooth = 1;

isfigure=0;
isoutput=1;
issyntest = 0;

%phvrange(1,:)=[3.55 4.15];
periods=2*pi./twloc;



lalim=[-11.2 -7.8];
lolim=[148.8 151.5];
gridsize=0.1;

%lalim=[30 50];
%lolim=[-125 -90];
%gridsize=0.3;

raydensetol=deg2km(gridsize)*3;
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
% for ie = 1
        for ip=1:length(periods)
%     for ip=6
        
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
        dt=-event(ie,ip).dt;
        if issyntest
            dt(:) = -event(ie,ip).ddist / 4;
        end
        if size(dt,1) == 1
            dt = dt';
        end
        para = polyfit(event(ie,ip).ddist(:),dt(:),1);
        avgv = abs(1./para(1));
        rays=event(ie,ip).ray;
        W = sparse(length(dt),length(dt));
        for i=1:length(dt)
            W(i,i)=1./event(ie,ip).fiterr(i);
            if issyntest
                W(i,i)=1;
            end
        end
        ind = find(W > maxerrweight);
        W(ind) = maxerrweight;
        ind = find(W < 1/fiterrtol);
        W(ind) = 0;
        for i=1:length(dt)
            %W(i,i)=W(i,i)*event(ie,ip).coherenum(i);
        end
        
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
        
        % Build the rotation matrix
        razi = azimuth(xi,yi,event(ie,ip).evla,event(ie,ip).evlo)+180;
        R = sparse(2*Nx*Ny,2*Nx*Ny);
        for i=1:Nx
            for j=1:Ny
                n=Ny*(i-1)+j;
                theta = razi(i,j);
                R(2*n-1,2*n-1) = cosd(theta);
                R(2*n-1,2*n) = sind(theta);
                R(2*n,2*n-1) = -sind(theta);
                R(2*n,2*n) = cosd(theta);
            end
        end
        
        % build dumping matrix for St
        dumpmatT = R(2:2:2*Nx*Ny,:);
        NR=norm(dumpmatT,1);
        NA=norm(W*mat,1);
        dumpweightT = Tdumpweight0*NA/NR;
        
        % build dumping matrix for SR
        dumpmatR = R(1:2:2*Nx*Ny-1,:);
        NR=norm(dumpmatR,1);
        NA=norm(W*mat,1);
        dumpweightR = Rdumpweight0*NA/NR;
        
        
        % build the smoothing operator
        smweight = smweight0;
        NR=norm(F,1);
        NA=norm(W*mat,1);
        smweight = smweight0*NA/NR;
        
        disp('start inverse');
        
        if isRsmooth
            A=[W*mat;smweight*F*R;dumpweightT*dumpmatT;dumpweightR*dumpmatR];
        else
            A=[W*mat;smweight*F;dumpweightT*dumpmatT;dumpweightR*dumpmatR];
        end
        rhs=[W*dt;zeros(size(F,1),1);zeros(size(dumpmatT,1),1);dumpweightR*ones(size(dumpmatR,1),1)./avgv];
        
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
            if stderr > dterrtol
                stderr = dterrtol;
            end
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
            
            % rescale dumping matrix for St
            NR=norm(dumpmatT,1);
            NA=norm(W*mat,1);
            dumpweightT = Tdumpweight0*NA/NR;
            
            % rescale dumping matrix for SR
            NR=norm(dumpmatR,1);
            NA=norm(W*mat,1);
            dumpweightR = Rdumpweight0*NA/NR;
            
            if isRsmooth
                A=[W*mat;smweight*F*R;dumpweightT*dumpmatT;dumpweightR*dumpmatR];
            else
                A=[W*mat;smweight*F;dumpweightT*dumpmatT;dumpweightR*dumpmatR];
            end
            rhs=[W*dt;zeros(size(F,1),1);zeros(size(dumpmatT,1),1);dumpweightR*ones(size(dumpmatR,1),1)./avgv];
            
            %            disp('start inverse');
            %            tic
            phaseg=(A'*A)\(A'*rhs);
            %            toc
            %            disp('Done');
            
        end
        
        % Calculate the kernel density
        %sumG=sum(abs(mat),1);
        ind=1:Nx*Ny;
        rayW = W;
        rayW(find(rayW>1))=1;
        raymat = rayW*mat;
        sumG(ind)=sum((raymat(:,2*ind).^2+raymat(:,2*ind-1).^2).^.5,1);
        clear raydense
        for i=1:Nx
            for j=1:Ny
                n=Ny*(i-1)+j;
                raydense(i,j)=sumG(n);
            end
        end
        
        %        disp(' Get rid of uncertainty area');
        for i=1:Nx
            for j=1:Ny
                n=Ny*(i-1)+j;
                if raydense(i,j) < raydensetol %&& ~issyntest
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
        if ~issyntest
            GV(ind) = NaN;
        end
        
        event_tomo(ie,ip).GV = full(GV);
        event_tomo(ie,ip).raydense = full(raydense);
        
        if isfigure
            figure(22)
            clf
            ax = worldmap(lalim, lolim);
            set(ax, 'Visible', 'off');
            surfacem(xi,yi,GV);
            load seiscmap
            colormap(seiscmap)
            drawpng
            colorbar
            avgphv = nanmean(raytomo(ip).GV(:));
            r = 0.2;
            caxis([avgphv*(1-r) avgphv*(1+r)])
            if issyntest
                caxis([3.5 4.5])
                for iray = 1:size(rays,1)
                     err = mat*phaseg - dt;
                    if abs(err(iray))>2
                        plotm([rays(iray,1),rays(iray,3)],[rays(iray,2) rays(iray,4)],'r')
                    else
                        plotm([rays(iray,1),rays(iray,3)],[rays(iray,2) rays(iray,4)],'k')
                    end
                end
            end
            plotm(event(ie,ip).evla,event(ie,ip).evlo,'rv','markersize',20);
            title('apparent phase V')
            
            
            figure(23)
            clf
            ax = worldmap(lalim, lolim);
            set(ax, 'Visible', 'off');
            if ~issyntest
                surfacem(xi,yi,GV-raytomo(ip).GV);
            else
                surfacem(xi,yi,GV-4);
            end
            load seiscmap
            colormap(seiscmap)
            drawpng
            colorbar
            title('Error to ray tomo');
            caxis([-0.5 0.5])
            disp(['Error to Ray tomo: ' ,num2str(nanmean(abs(GV(:)-raytomo(ip).GV(:))))]);
            
            Sr = R*phaseg;
            for i=1:Nx
                for j=1:Ny
                    n=Ny*(i-1)+j;
                    GVr(i,j)= Sr(2*n-1);
                    GVt(i,j)= Sr(2*n);
                end
            end
            figure(24)
            clf
            ax = worldmap(lalim, lolim);
            set(ax, 'Visible', 'off');
            surfacem(xi,yi,GVr);
            load seiscmap
            colormap(seiscmap)
            drawpng
            colorbar
            title('Sr');
            figure(25)
            clf
            ax = worldmap(lalim, lolim);
            set(ax, 'Visible', 'off');
            surfacem(xi,yi,GVt);
            load seiscmap
            colormap(seiscmap)
            drawpng
            colorbar
            title('St');
        end
        
    end	% end of period loop
end % end of event loop
if isoutput
    save('event_tomo.mat','event_tomo','xnode','ynode','periods');
end