function [G1]=DFT_2D(X_pos,Y_pos,Kx,Ky,S,delta_k)

    G2=zeros(length(X_pos),length(Y_pos));%define the Fourier series
    for kk=1:length(Y_pos)
        for jj=1:length(X_pos)
            for ii=1:length(Kx)
                for mm=1:length(Ky)

                    chi=Kx(ii)*X_pos(jj)+Ky(mm)*Y_pos(kk);
                    G2(jj,kk)=exp(-1i*2*pi*chi)*S(ii,mm)+G2(jj,kk);
                end
            end
        end
    end
        G1=(delta_k/(2*pi))^2.*G2/10^6;
end