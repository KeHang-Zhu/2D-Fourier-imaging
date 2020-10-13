function [B,phi,phase,imaging]=Mag_senspr_2D(N,B0,HW,pos_NVx,pos_NVy,int,n_spin,gamma,T2,I)
%     N=Size/pixel;
    Conv=1000;
    x=(1:N);
    x=x/Conv;%change into um
 if I==1
    %% set the lorenzian magnetic field profile (giant peak)
    B=zeros(N,N);
    for ii=1:N
        for jj=1:N
             B(ii,jj)=B0/(1+((ii-N/2)*int/HW)^2)/(1+((jj-N/2)*int/HW)^2);
        end
    end
    %     figure
    %     mesh(x,x,B);
    %     view(2)
    %     xlabel('x/um');
    %     ylabel('y/um');
    %     title('Field ditribution');
    %     colorbar
    %
 elseif I==2
     %% set the dipolar magnetic field profile  
     B=0;zzz=30*int;
     for num=1:length(pts)
         xc = pts(num,1)*int;
         yc = pts(num,2)*int;
         
         for ii=1:NN*N
             for jj=1:NN*N
                 %val(ii,jj)=B0/(1+((ii-xc)*int/HW)^2)/(1+((jj-yc)*int/HW)^2);
                 r=sqrt(zzz^2+(ii*int-xc)^2+(jj*int-yc)^2);
                 costheta=zzz/r;
                 val(ii,jj)=(3*costheta^2-1)*u/r^3;
             end
         end
         
         B=B+val;
         
     end
     BB=B(5*N+1:6*N,5*N+1:6*N);
     figure
     surf(x,y,BB')
     
     shading interp
     view(2)
     colorbar
     xlabel('x/um');
     ylabel('y/um');
     title('magnetic field ditribution');
     for ii=1:N
         for jj=1:N
             BB(ii,jj)=exp(1i*BB(ii,jj));
         end
     end
 end
 %% plotiong
    figure
    hold on
    for ii=1:1:n_spin
         for jj=1:1:n_spin
            plot3(pos_NVx(ii,jj)/Conv,pos_NVy(ii,jj)/Conv,1,'r.')
         end
    end
    hold off
    xlabel('x/um');
    ylabel('y/um');
    title('sensor ditribution');
    
%%  Phase ditribution
    phi=zeros(N);
    for ii=1:N
        for jj=1:N
            %                              phi((ii-1)*n_read+jj)=exp(1i*2*pi*gamma*T2*B(num));
            phi(ii,jj)=gamma*T2*B(ii,jj);% initial phase
        end
    end
    
    figure
    surf(x,x,phi);
    view(2);
    colorbar
    xlabel('x/um');
    ylabel('y/um');
    title('Phase ditribution');
%% point loss
    phase=zeros(length(phi(1,:)));
    imaging=zeros(length(phi(1,:)));
    AAA=pos_NVx(:); BBB=pos_NVy(:);
    for ii=1:length(AAA)
            phase(AAA(ii),BBB(ii))=phi(AAA(ii),BBB(ii));
            imaging(AAA(ii),BBB(ii))=exp(1i*2*pi*phase(AAA(ii),BBB(ii)));
    end
    figure
    surf(x,x,phase)
    grid off
    %shading interp
    view(2)
    hold off
    xlabel('x/um');
    ylabel('y/um');
    title('field ditribution');
end