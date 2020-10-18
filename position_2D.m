 function [pos_NVx,pos_NVy]=position_2D(N,n_spin,distance,disturbed)
% function [pts]=position_2D(N,n_spin,distance,mode,disturbed)
% mode stands for what kind of 2D lattice 
% disturbed stands for whether the lattice is perfect or not

      pos_NVx=zeros(n_spin,n_spin);
      pos_NVy=zeros(n_spin,n_spin);
% if I==1
%     ranx=rand(n_spin,n_spin);%Add perturbation to the NV lattice
%     rany=rand(n_spin,n_spin);%Add perturbation to the NV lattice
%     for jj=1:1:n_spin
%         for ii=1:1:n_spin %NV's position relative to the magnetic field
%             pos_NVx(ii,jj)=ceil(1+distance*(ii-1)+sqrt(distance)*2*ranx(ii,jj));
%             pos_NVy(ii,jj)=ceil(1+distance*(jj-1)+sqrt(distance)*2*rany(ii,jj));
%             if  pos_NVy(ii,jj)>N
%                 pos_NVy(ii,jj)=N;
%             end
%             
%             if  pos_NVx(ii,jj)>N
%                 pos_NVx(ii,jj)=N;
%             end
%         end
%     end
%     
% if mode==0
    if disturbed==0
        ranx=zeros(n_spin,n_spin);% ideal case
        rany=zeros(n_spin,n_spin);% ideal case
    elseif disturbed==1
        ranx=rand(n_spin,n_spin);%Add perturbation to the NV lattice 
        rany=rand(n_spin,n_spin);%Add perturbation to the NV lattice 
    end
        
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

% spacing=distance;
%     
% if mode==1
%     
%      pts =squareGrid([0 0 N N], [N/2 N/2], [spacing,spacing]);
%      
%      if disturbed==0
%          ranx=zeros(length(pts),1);% ideal case
%          rany=zeros(length(pts),1);% ideal case
%      elseif disturbed==1
%          ranx=rand(length(pts),1)*sqrt(distance);%Add perturbation to the NV lattice
%          rany=rand(length(pts),1)*sqrt(distance);%Add perturbation to the NV lattice
%      end
%      
%      pts(:,1)=pts(:,1)+ranx;
%      pts(:,2)=pts(:,2)+rany;
%      
%      for ii=1:length(pts)
%          if  pts(ii,1)>N
%              pts(ii,1)=N;
%          end
%          if  pts(ii,2)>N
%              pts(ii,2)=N;
%          end
%          if  pts(ii,1)<1
%              pts(ii,1)=1;
%          end
%          if  pts(ii,2)<1
%              pts(ii,2)=1;
%          end
%      end
% 
%     figure
%     scatter(pts(:,1),pts(:,2))
%     xlabel('x/um');
%     ylabel('y/um');
%     title('rectangular superlattice');
%      
% 
% elseif mode==2
%     %% set the triangular-lattice
%     %spacing=distance;
%     %NN=11;  %larger size to reduce the marginal effect of finite sample size
%     
%     %  pts =squareGrid([0 0 3*N 3*N], [3*N/2 3*N/2], [spacing,spacing]);
%     pts =triangleGrid([0 0 N N], [N/2 N/2], spacing);
%     
%     if disturbed==0
%         ranx=zeros(1,length(pts));% ideal case
%         rany=zeros(1,length(pts));% ideal case
%     elseif disturbed==1
%         ranx=rand(1,length(pts))*sqrt(distance);%Add perturbation to the NV lattice
%         rany=rand(1,length(pts))*sqrt(distance);%Add perturbation to the NV lattice
%     end
%     
%      pts(:,1)=pts(:,1)+ranx;
%      pts(:,2)=pts(:,2)+rany;
%      
%      for ii=1:length(pts)
%          if  pts(ii,1)>N
%              pts(ii,1)=N;
%          end
%          if  pts(ii,2)>N
%              pts(ii,2)=N;
%          end
%          if  pts(ii,1)<1
%              pts(ii,1)=1;
%          end
%          if  pts(ii,2)<1
%              pts(ii,2)=1;
%          end
%      end
%     %% plot the NV lattice
%     figure
%     scatter(pts(:,1),pts(:,2))
%     xlabel('x/um');
%     ylabel('y/um');
%     title('triangular superlattice');
%     
%     %% plot the NV lattice
% %     latt=zeros(NN*N+1);
% %     for ii=1:length(pts)
% %         latt(round(pts(ii,1))+1,round(pts(ii,2))+1)=1;
% %     end
% %     latt100=latt(5*N+1:6*N,5*N+1:6*N);
% %     figure
% %     mesh(1:1:100,1:1:100,latt100')
%     
% elseif mode==3
%       %% set the hex-lattice
%        %NN=11;  %larger size to reduce the marginal effect of finite sample size
%        pts =hexagonalGrid([0 0 N N], [N/2 N/2], spacing);
%        
%        if disturbed==0
%            ranx=zeros(1,length(pts));% ideal case
%            rany=zeros(1,length(pts));% ideal case
%        elseif disturbed==1
%            ranx=rand(1,length(pts))*sqrt(distance);%Add perturbation to the NV lattice
%            rany=rand(1,length(pts))*sqrt(distance);%Add perturbation to the NV lattice
%        end
%        
%        pts(1,:)=pts(1,:)+ranx;
%        pts(2,:)=pts(2,:)+rany;
%        
%        for ii=1:length(pts)
%            if  pts(ii)>N
%                pts(ii)=N;
%            end
%        end
%      %%
%     figure
%     scatter(pts(:,1),pts(:,2))
%     xlabel('x/um');
%     ylabel('y/um');
%     title('hexigonal superlattice');
     %%
%     latt=zeros(NN*N+1);
%     for ii=1:length(pts)
%         latt(round(pts(ii,1))+1,round(pts(ii,2))+1)=1;
%     end
%     latt100=latt(5*N+1:6*N,5*N+1:6*N);
%     figure
%     mesh(1:1:100,1:1:100,latt100')
% end
end