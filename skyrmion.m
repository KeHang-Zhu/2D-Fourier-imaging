function [BB]=skyrmion(N,z,r)

xc=N/2;
yc=N/2;
pts=zeros(4,1);
pts(1,1)=xc;
pts(2,1)=yc;
pts(3,1)=-pi;
pts(4,1)=0;
ub=9.274*10^(-24);
u=ub*10^(-4);%magnetism for spin 1/2 particle

number=6;
radius=r;
int=10^(-9);
zzz=z*int;

for ii=1:4
    theta=2*pi/number/2^(ii-1)/2:2*pi/number/2^(ii-1):2*pi;
    for jj=1:length(theta)
        the=theta(jj);
        pts(1,end+1)=cos(the)*radius*ii+xc;
        pts(2,end)=sin(the)*radius*ii+yc;
        pts(3,end)=-pi+pi/4*ii;
        pts(4,end)=the;
    end
end
% figure
% scatter(pts(1,:),pts(2,:))
% set(gca,'XLim',[0,N],'YLim',[0,N])
% xlabel('x/um');
% ylabel('y/um');
% title('Skyrmion spin configuration');

BB=0;
val=zeros(N,N);

for nu=1:length(pts(1,:))
    xa = pts(1,nu)*int;
    ya = pts(2,nu)*int;
    mz=cos(pts(3,nu));
    mx=sin(pts(3,nu))*cos(pts(4,nu));
    my=sin(pts(3,nu))*sin(pts(4,nu));
    for ii=1:N
        for jj=1:N
            
            r=sqrt(zzz^2+(ii*int-xa)^2+(jj*int-ya)^2);
            rx=(ii*int-xa)/r;
            ry=(jj*int-ya)/r;
            rz=zzz/r;
            
            costheta=zzz/r;
            rdm=rx*mx+ry*my+rz*mz;
            val(ii,jj)=(-mz+3*costheta*rdm)*u/r^3;
            
        end
    end
    
    BB=BB+val;
    
end

x=1/1000:1/1000:N/1000;

figure
mesh(x,x,BB);
view(2);
colorbar
xlabel('x/um');
ylabel('y/um');
title('Phase ditribution');
end