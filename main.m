%% parameter
a0=1*10^(-8);%NV distance
aa0=a0*10^(6);%convert the average distance of NV into um
HW=8*10^(-8);%half width of the magnetic field distribution
B0=0.2;%amplitude of the magnetic field Gauss ,therefore phi_s(x)=2pi*gamma*T2*B_max won't exceeds 2pi, avoids the dillemma
gamma=2.803;%gyromagnetic ratio
int=10^(-9);%pixel size 1 nm
Size=10*10^(-7);%sample size 100nm
n_read=Size/a0;%readout 11 NV at once
N=100;%number of points
Bg=20;%max gradient 
T2=0.3;%T2* 0.3us
C=1;%contrast 1
distance=ceil(a0/int);
n_spin=ceil(N/distance);%number of spins
space=n_spin/n_read;
delta_x=a0*10^6;%distance between NV centers in scale of um
delta_k=1/Size/10^(6);%field of view
delta_y=a0*10^6;%distance between NV centers in scale of um
conv2um=1000;
%% %NV's position 
[pos_NVx,pos_NVy]=position_2D(N,n_spin,distance,2);%0 represent regular rectangular lattice, 
                                                    %1 represents distrubed rectangular lattice
                                                    %2 represents regular triangular lattice, 
                                                    %3 represents regular hex lattice                                                    
%% scanning part    

spacing=50; %real space scanning 50 nm
N_px=3;%sampling point in x direction
N_py=3;%sampling point in y direction
Recon=zores(N_px,N_py,N*N);
scanning_x=zores(N_px,N_py,N*N);
scanning_y=zores(N_px,N_py,N*N);
Cont=zeros(length(N_p),length(Kmax));

%% set the field at the sensor location
[MB]=Mag_senspr_2D(N,B0,HW,pos_NVx,pos_NVy,int,n_spin,gamma,T2,0);
                                                    %0 represent single broad peak, 
                                                    %1 represent dipolar field generate by ferromagnetical aligned spins

for jj=1:length(N_px)
 for mm=1:length(N_py)

B=MB((jj-1)*spacing:(jj-1)*spacing+N,(mm-1)*spacing:(mm-1)*spacing+N); 

%% Ramsey sequence
Kxmax=N;
Kymax=Kxmax;
Kx=1:1:Kxmax ; %x wave vector
Ky=1:1:Kymax ; %y wave vector
Fs = 1/Kx(1); %field of view
[S,phi]=ksample(N,pos_NVx,pos_NVy,n_spin,B,T2,gamma,Kxmax,Kymax,0,[]);
% Im=imag(S);
% Re=real(S);
   
%% inverse Fourier transformation
[G2]=invFourier_2D(S,N,N);
ABS=abs(G2);
phase=angle(G2);
 
pixel=1/1000;%pixel size=1000/diatance/NNN nm
pos_x=pixel:pixel:N/1000;
pos_y=pixel:pixel:N/1000;
figure
mesh(pos_x,pos_y,ABS);
view(2)
colorbar
xlabel('x/um');
ylabel('y/um');
title('reconstructed NV center ditribution'); 
 
figure
mesh(pos_x,pos_y,phase);
view(2)
colorbar
xlabel('x/um');
ylabel('y/um');
title('reconstructed field ditribution'); 
 
 %% track the B distribution after locating the NV distribution
 %ploting_2D(N,n_spin,X,Y,ABS);
[pos_Xr,pos_Yr,X_r,Y_r,phase_re]=locating_2D(ABS,phase,n_spin,pos_x,pos_y); 

%% collect scanning points
Recon(jj,mm,1:size(phase_re))=phase_re(:);
X_r=X_r+(jj-1)*spacing/con2um;
Y_r=Y_r+(jj-1)*spacing/con2um;
scanning_x(jj,mm,1:length(X_r))=X_r(:);
scanning_y(jj,mm,1:length(Y_r))=Y_r(:);
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
import2R=[XX./con2um,YY./con2um,val];

figure
mesh(import2R(:,1),import2R(:,2),import2R(:,3));
view(2)
colorbar
xlabel('x/um');
ylabel('y/um');
title('scanned reconstruction field'); 

%% curve fitting for getting the info of dipolar kernal
% find the suitable initial values for curve fitting
B0=max(val);



 