function [G2]=invFourier_2D(S,Kxmax,Kymax)

    G2=ifft2(S,Kxmax,Kymax);  
end