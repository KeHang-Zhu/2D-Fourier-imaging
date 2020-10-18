function [X_r,Y_r,X_r1,Y_r1,phase]=locating_2D(ABS,PHA,n_spin,x,y)
    conv=1000;
    benchmark=0.7;
    [pks,X_r,Y_r] = findpeak_2D(ABS);
    %[locs]=FastPeakFind(ABS);
%     X_r=zeros(1,length(locs)/2);
%     Y_r=zeros(1,length(locs)/2);
%     pks=zeros(1,length(locs)/2);
%     for ii=1:length(locs)/2
%         X_r(ii)=locs(ii);
%         Y_r(ii)=locs(2*ii);
%         pks(ii)=ABS(X_r(ii),Y_r(ii));
%     end

    for ii=length(pks):-1:1
        if pks(ii)<=benchmark
            X_r(ii)=[];
            Y_r(ii)=[];
        end
    end
       
    X_r1=X_r./conv;% convert the scale to um
    Y_r1=Y_r./conv;
 %% plotting the located NV centers
    figure
    hold on
    for ii=1:1:length(X_r)
         
            plot3(X_r1(ii),Y_r1(ii),1,'r.');
        
    end
    hold off
    xlabel('x/um');
    ylabel('y/um');
    title('reconstructed sensor ditribution');
    
%% plotting reconstructed field ditribution after locating the NV center
    phase=zeros(length(PHA));
    for ii=1:length(X_r)
            phase(X_r(ii),Y_r(ii))=PHA(X_r(ii),Y_r(ii));
    end
    figure
    mesh(x,y,phase)
    grid off
    shading interp
    view(2)
    hold off
    xlabel('x/um');
    ylabel('y/um');
    title('reconstructed field ditribution after locating the NV center');

    
end