function [pos_NVx,pos_NVy]=position_2D(N,n_spin,distance,I)
      pos_NVx=zeros(n_spin,n_spin);
       pos_NVy=zeros(n_spin,n_spin);
%       ranx=rand(n_spin,n_spin);%Add perturbation to the NV lattice 
%       rany=rand(n_spin,n_spin);%Add perturbation to the NV lattice 
%       ranx=zeros(n_spin,n_spin);% ideal case  
%       rany=zeros(n_spin,n_spin);% ideal case
if I==1
      ranx=rand(n_spin,n_spin);%Add perturbation to the NV lattice 
      rany=rand(n_spin,n_spin);%Add perturbation to the NV lattice 
    for jj=1:1:n_spin
        for ii=1:1:n_spin %NV's position relative to the magnetic field
            pos_NVx(ii,jj)=ceil(1+distance*(ii-1)+sqrt(distance)*2*ranx(ii,jj));
            pos_NVy(ii,jj)=ceil(1+distance*(jj-1)+sqrt(distance)*2*rany(ii,jj));
            if  pos_NVy(ii,jj)>N
                pos_NVy(ii,jj)=N;
            end
            
            if  pos_NVx(ii,jj)>N
                pos_NVx(ii,jj)=N;
            end            
        end
    end
    
elseif I==0
    ranx=zeros(n_spin,n_spin);% ideal case
    rany=zeros(n_spin,n_spin);% ideal case
     for jj=1:1:n_spin
        for ii=1:1:n_spin %NV's position relative to the magnetic field
            pos_NVx(ii,jj)=ceil(1+distance*(ii-1)+sqrt(distance)*2*ranx(ii,jj));
            pos_NVy(ii,jj)=ceil(1+distance*(jj-1)+sqrt(distance)*2*rany(ii,jj));
            if  pos_NVy(ii,jj)>N
                pos_NVy(ii,jj)=N;
            end
            
            if  pos_NVx(ii,jj)>N
                pos_NVx(ii,jj)=N;
            end   
        end
    end


elseif I==2
    %% set the triangular-lattice
    spacing=distance;
    NN=11;  %larger size to reduce the marginal effect of finite sample size
    
    %  pts =squareGrid([0 0 3*N 3*N], [3*N/2 3*N/2], [spacing,spacing]);
    pts =triangleGrid([0 0 NN*N NN*N], [NN*N/2 NN*N/2], spacing);
    %%
    figure
    scatter(pts(:,1),pts(:,2))
    xlabel('x/um');
    ylabel('y/um');
    title('triangular superlattice');
    %% rotation
    % angl=pi/5;
    % for ii=1:length(pts(:,1))
    %     pts(ii,1)=pts(ii,1)*cos(angl)-pts(ii,2)*sin(angl);
    %     pts(ii,2)=pts(ii,2)*cos(angl)+pts(ii,1)*sin(angl);
    % end
    %%
    latt=zeros(NN*N+1);
    for ii=1:length(pts)
        latt(round(pts(ii,1))+1,round(pts(ii,2))+1)=1;
    end
    latt100=latt(5*N+1:6*N,5*N+1:6*N);
    figure
    mesh(1:1:100,1:1:100,latt100')
elseif I==3
      %% set the hex-lattice
      spacing=distance;
       NN=11;  %larger size to reduce the marginal effect of finite sample size
       pts =hexagonalGrid([0 0 NN*N NN*N], [NN*N/2 NN*N/2], spacing);
     %%
    figure
    scatter(pts(:,1),pts(:,2))
    xlabel('x/um');
    ylabel('y/um');
    title('hexigonal superlattice');
     %%
    latt=zeros(NN*N+1);
    for ii=1:length(pts)
        latt(round(pts(ii,1))+1,round(pts(ii,2))+1)=1;
    end
    latt100=latt(5*N+1:6*N,5*N+1:6*N);
    figure
    mesh(1:1:100,1:1:100,latt100')
end
end