function []=ploting_2D(N,n_spin,X,Y,phase,pos_Xr,pos_Yr,NNN,pixel)
    PHA=zeros(n_spin*NNN);
    posx=round(pos_Xr./pixel);
    posy=round(pos_Yr./pixel);
    for ii=length(posx)
        for jj=1:length(posx)
            n1=posx(ii,jj);
            n2=posy(ii,jj);
            PHA(n1,n2)=phase(ii,jj);
        end
    end
    figure
    surf(X,Y,PHA);
    view(2);
    xlabel('x/um');
    ylabel('y/um');
    title('phi distribution');
   colorbar
    
end