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
HW=8*10^(-8);%half width of the magnetic field distribution
B0=0.2;%amplitude of the magnetic field Gauss ,therefore phi_s(x)=2pi*gamma*T2*B_max won't exceeds 2pi, avoids the dillemma
Size=10*10^(-7);%sample size 100nm
period=20;  %spatial period of magnetic field
delta_k=1/Size/10^(6);%field of view

%% lattice playground
disturbed=0;      %whether the NV grid is ideal or random

%underlying super-lattice structure
mode=1;                 %1 represents rectangular lattice
                        %2 represents regular triangular lattice, 
                        %3 represents regular hex lattice  
field=1;             %1 represent single broad lorenzian peak (for test), 
                     %2 represent dipolar field generate by ferromagnetical aligned spins        
                                                 
%% scanning part    

spacing=200; %real space scanning 50 nm
N_px=1;%scanning points in x direction
N_py=1;%scanning points in y direction
Recon=zeros(N_px,N_py,N*N);
scanning_x=zeros(N_px,N_py,N*N);
scanning_y=zeros(N_px,N_py,N*N);

%% %NV's position 
 [pos_NVx,pos_NVy]=position_2D(N,n_spin,distance,disturbed); 
 
%% set the field at the sensor location (with wider region)
[MB]=Mag_senspr_2D(N,B0,HW,pos_NVx,pos_NVy,int,n_spin,period,height,gamma,T2,field,mode);   

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
     % Im=imag(S);
     % Re=real(S);   
     
     %% inverse Fourier transformation
     [G2]=invFourier_2D(S,N,N);
     ABS=abs(G2);
     phase=angle(G2);
     phase=phase./(2*pi*T2*gamma);
   
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
     %ploting_2D(N,n_spin,X,Y,ABS);
     [pos_Xr,pos_Yr,X_r,Y_r,phase_re]=locating_2D(ABS,phase,n_spin,pos_x,pos_y);
     
     %% collect scanning points
     Recon(jj,mm,1:size(phase_re(:)))=phase_re(:);
     X_r=X_r+(jj-1)*spacing/con2um;
     Y_r=Y_r+(mm-1)*spacing/con2um;
     scanning_x(jj,mm,1:length(X_r(:)))=X_r(:);
     scanning_y(jj,mm,1:length(Y_r(:)))=Y_r(:);
  end
end

%% reshape real space scanning data
val=Recon(:);
XX=scanning_x(:);
YY=scanning_y(:);
% XX(find(XX==0))=[];%take out the zero 
% YY(find(YY==0))=[];
% [XX,index]=unique(XX); !take out the repeated points
% yy=zeros(length(XX),1);
% for ii=1:length(XX)
%     yy(ii)=YY(index(ii));
% end

%% generate the data form convenient for R language processing
import2R=[XX,YY,val];

% figure
% plot3(import2R(:,1),import2R(:,2),import2R(:,3));
% view(2)
% colorbar
% xlabel('x/um');
% ylabel('y/um');
% title('scanned reconstruction field'); 

%% curve fitting for getting the info of dipolar kernal
% find the suitable initial values for curve fitting
[B0,index_c]=max(val);
xc=XX(index_c);
yc=YY(index_c);

%% cluster all the points into one unit cell
for ii=1:length(XX)
    deltax=XX(ii)-xc;
    deltay=YY(ii)-yc;
    re_x=rem(deltax,period);
    re_y=rem(deltay,period);
end
  import2R_cluster=[re_x,re_y,val];  




 