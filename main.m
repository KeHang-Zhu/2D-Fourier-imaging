%% constant
gamma=2.803;%gyromagnetic ratio
ub=9.274*10^(-24);
u=ub;%magnetism for spin 1/2 particle
conv2um=1000;
int=10^(-9);%pixel size 1 nm
T2=0.3;%T2* 0.3us
N=100;%number of points

%% parameter NV related
a0=1*10^(-8);%NV distance
aa0=a0*10^(6);%convert the average distance of NV into um
distance=ceil(a0/int);
n_spin=ceil(N/distance);%number of spins
delta_x=a0*10^6;%distance between NV centers in scale of um
delta_y=a0*10^6;%distance between NV centers in scale of um
C=1;%contrast 1

%% parameter material related
height=30;%the distance between NV center and material
HW=2*10^(-8);%half width of the magnetic field distribution
B0=0.2;%amplitude of the magnetic field Gauss ,therefore phi_s(x)=2pi*gamma*T2*B_max won't exceeds 2pi, avoids the dillemma
Size=10*10^(-7);%sample size 100nm
period=20;  %spatial period of magnetic field
delta_k=1/Size/10^(6);%field of view

%% lattice playground
disturbed=1;      %whether the NV grid is ideal or random

%underlying super-lattice structure
mode=1;                 %1 represents rectangular lattice
                        %2 represents regular triangular lattice, 
                        %3 represents regular hex lattice  
field=2;             %1 represent single broad lorenzian peak (for test), 
                     %2 represent dipolar field generate by ferromagnetical aligned spins        
                                                 
%% scanning part    
spacing=50; %real space scanning 50 nm
N_px=2;%scanning points in x direction
N_py=2;%scanning points in y direction
Recon=zeros(N_px,N_py,n_spin*n_spin);
scanning_x=zeros(N_px,N_py,n_spin*n_spin);
scanning_y=zeros(N_px,N_py,n_spin*n_spin);

%% %NV's position 
 [pos_NVx,pos_NVy]=position_2D(N,n_spin,distance,disturbed); 
 
%% set the field at the sensor location (with wider region)
 mul=spacing*(max(N_px,N_py)-1)+N;
 [MB]=Mag_senspr_2D(N,mul,B0,HW,pos_NVx,pos_NVy,int,n_spin,period,height,gamma,T2,field,mode); 
 
 pos_x=pixel:pixel:N/1000;
 pos_y=pixel:pixel:N/1000;

for jj=1:N_px
  for mm=1:N_py

     B=MB((jj-1)*spacing+1:(jj-1)*spacing+N,(mm-1)*spacing+1:(mm-1)*spacing+N); % the B field at the moved region 
     
%      pixel=1/1000;%pixel size=1000/diatance/NNN nm
%      pos_x=pixel:pixel:N/1000;
%      pos_y=pixel:pixel:N/1000;
%      figure
%      mesh(pos_x,pos_y,B);
%      view(2)
%      colorbar
%      xlabel('x/um');
%      ylabel('y/um');
%      title('input magnetic field ditribution');
     
     %% Ramsey sequence
     Kxmax=N;
     Kymax=Kxmax;
     Kx=1:1:Kxmax ; %x wave vector
     Ky=1:1:Kymax ; %y wave vector
     Fs = 1/Kx(1); %field of view
     [S,phi]=ksample(N,pos_NVx,pos_NVy,n_spin,B,T2,gamma,Kxmax,Kymax,0,[]); %perform Ramsey sequence measurement

      %%  plot the K space spectrum  
%     ABS=abs(S);
%     Re=real(S);
%     deltak=conv2um/N;
%     fx=deltak:deltak:1000;
%     fy=fx;
%     figure
%     mesh(fx,fy,ABS);
%     view(2)
%     colorbar
%     xlabel('kx/um^{-1}');
%     ylabel('ky/um^{-1}');
%     title('F(kx,ky)');

    %% compressed sensing (by randomly tossing S points)
%     Tosspt=floor(N/5*4);
%     Tosspt_x=sort(randperm(N,Tosspt));
%     Tosspt_y=sort(randperm(N,Tosspt)); 
%     for ii=Tosspt:-1:1
%         for ll=Tosspt:-1:1
%             x_t=Tosspt_x(ii);
%             y_t=Tosspt_y(ll);
%             S(x_t,y_t)=0;
%         end
%     end
    
     %%  plot the K space spectrum of compressed sensing
%     ABS=abs(S);
%     deltak=conv2um/N;
%     fx=deltak:deltak:1000;
%     fy=fx;
%     figure
%     mesh(fx,fy,ABS);
%     view(2)
%     colorbar
%     xlabel('kx/um^{-1}');
%     ylabel('ky/um^{-1}');
%     title('compressed sensing');
    
     %% inverse Fourier transformation
     [G2]=invFourier_2D(S,N,N);
     ABS=abs(G2);
     phase=angle(G2);
     phase=phase./(2*pi*T2*gamma);
     
     %% plot the original inv Fourier pattern
%      pixel=1/1000;%pixel size=1000/diatance/NNN nm
%      pos_x=pixel:pixel:N/1000;
%      pos_y=pixel:pixel:N/1000;
%      figure
%      mesh(pos_x,pos_y,ABS);
%      view(2)
%      colorbar
%      xlabel('x/um');
%      ylabel('y/um');
%      title('reconstructed NV center ditribution');
%      
%      figure
%      mesh(pos_x,pos_y,phase);
%      view(2)
%      colorbar
%      xlabel('x/um');
%      ylabel('y/um');
%      title('reconstructed field ditribution');
     
     %% track the B distribution after locating the NV distribution
     
     [pos_Xr,pos_Yr,X_r,Y_r,phase_re,non_pha]=locating_2D(ABS,phase,n_spin,pos_x,pos_y);
     
     %% collect scanning points
     Recon(jj,mm,1:length(X_r))=non_pha;
     X_r=X_r+(jj-1)*spacing/con2um;
     Y_r=Y_r+(mm-1)*spacing/con2um;
     scanning_x(jj,mm,1:length(X_r))=X_r;
     scanning_y(jj,mm,1:length(Y_r))=Y_r;
  end
end

%% reshape real space scanning data
val=Recon(:);
XX=scanning_x(:);
YY=scanning_y(:);

% [XX,index]=unique(XX); %take out the repeated points
% yy=zeros(length(XX),1);
% val2=zeros(length(XX),1);
% for ii=1:length(XX)
%     yy(ii)=YY(index(ii));
%     val2(ii)=val(index(ii));
% end

%% generate the data form convenient for R language processing
import2R=[XX,YY,val];

% figure
% plot3(import2R(:,1),import2R(:,2),import2R(:,3),'.');
% xlabel('x/um');
% ylabel('y/um');
% title('scanned reconstruction field'); 

if (field==2)
    %% curve fitting for getting the info of dipolar kernal
    % find the suitable initial values for curve fitting
    [B0,index_c]=max(val);
    xc=XX(index_c);
    yc=YY(index_c);
    
    %% cluster all the points into one unit cell
    re_x=XX;re_y=YY;
    for ii=1:length(XX)
        deltax=XX(ii)-xc;
        deltay=YY(ii)-yc;
        re_x(ii)=rem(deltax,period/conv2um);
        re_y(ii)=rem(deltay,period/conv2um);
        if (re_x(ii)<-period/conv2um/2)
            re_x(ii)=re_x(ii)+period/conv2um;
        end
        if (re_x(ii)>period/conv2um/2)
            re_x(ii)=re_x(ii)-period/conv2um;
        end
        if (re_y(ii)<-period/conv2um/2)
            re_y(ii)=re_y(ii)+period/conv2um;
        end
        if (re_y(ii)>period/conv2um/2)
            re_y(ii)=re_y(ii)-period/conv2um;
        end
    end
    import2R_cluster=[re_x,re_y,val];
    
    
    figure
    plot3(import2R_cluster(:,1),import2R_cluster(:,2),import2R_cluster(:,3),'.');
    set(gca,'XLim',[-0.05 0.05],'YLim',[-0.05 0.05]);
    xlabel('x/um');
    ylabel('y/um');
    title('clusterred reconstruction field');
end
 