function [B,phi]=Mag_senspr_2D(N,mul,B0,HW,pos_NVx,pos_NVy,int,n_spin,distance,height,gamma,T2,I,J)
%     N=Size/pixel;
    Conv=1000;
    ub=9.274*10^(-24);
    u=ub;%magnetism for spin 1/2 particle
    pos_NVx=pos_NVx*int;
    pos_NVy=pos_NVy*int;
    x=(1:mul);
    x=x/Conv;%change into um
    
    
 if I==1
    %% set the lorenzian magnetic field profile (giant peak)
    B=zeros(mul,mul);
    for ii=1:mul
        for jj=1:mul
             B(ii,jj)=B0/(1+((ii-mul/2)*int/HW)^2)/(1+((jj-mul/2)*int/HW)^2);
        end
    end
    
%% plot the field
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
     BB=0;zzz=height*int;
     NN=11; %multiple of points to overcome the finite size effect
     NN2=NN/2;
     add=(mul)/N;
     add2=add/2;
     val=zeros(N*NN);
     
     spacing=distance;
     
     if J==1
         pts =squareGrid([0 0 NN*N NN*N], [NN*N/2 NN*N/2], [spacing,spacing]);  
     elseif J==2
         pts =triangleGrid([0 0 NN*N NN*N], [NN*N/2 NN*N/2], spacing);
     elseif J==3
         pts =hexagonalGrid([0 0 NN*N NN*N], [NN*N/2 NN*N/2], spacing);
     end
     
     %% plot
     figure
     scatter(pts(:,1),pts(:,2))
     xlabel('x/um');
     ylabel('y/um');
     title('triangular superlattice');
     
     for num=1:length(pts(:,1))
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
         
         BB=BB+val;
         
     end
     B=BB((NN2-add2)*N+1:(NN2+add2)*N,(NN2-add2)*N+1:(NN2+add2)*N);% we take the middle points to overcome the finite size effect

    %% plot the field
%      figure
%      surf(x,y,BB')
%      
%      shading interp
%      view(2)
%      colorbar
%      xlabel('x/um');
%      ylabel('y/um');
%      title('magnetic field ditribution');
%      for ii=1:N
%          for jj=1:N
%              BB(ii,jj)=exp(1i*BB(ii,jj));
%          end
%      end
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
    
    %% Phase ditribution
    phi=zeros(mul);
    for ii=1:mul
        for jj=1:mul
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
%     phase=zeros(length(phi(1,:)));
%     imaging=zeros(length(phi(1,:)));
%     AAA=pos_NVx(:); BBB=pos_NVy(:);
%     for ii=1:length(AAA)
%             phase(AAA(ii),BBB(ii))=phi(AAA(ii),BBB(ii));
%             imaging(AAA(ii),BBB(ii))=exp(1i*2*pi*phase(AAA(ii),BBB(ii)));
%     end
%     figure
%     surf(x,x,phase)
%     grid off
%     %shading interp
%     view(2)
%     hold off
%     xlabel('x/um');
%     ylabel('y/um');
%     title('field ditribution');
end