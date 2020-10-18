function [S_loss,phase]=ksample(N,pos_NVx,pos_NVy,n_spin,B,T2,gamma,Kxmax,Kymax,I,photon)%N is the multiple of 
    phi=zeros(n_spin,n_spin);
    pha=zeros(N,N);
    phase=zeros(N,N);
    for ii=1:n_spin
        for jj=1:n_spin
            num1=pos_NVx(ii,jj);
            num2=pos_NVy(ii,jj);
            phi(ii,jj)=gamma*T2*B(num1,num2);% initial phase
            pha(num1,num2)=exp(1i*2*pi*phi(ii,jj));
            phase(num1,num2)=2*pi*phi(ii,jj);
        end
    end

    S=fft2(pha,Kxmax,Kymax);
    
    %% add noise (photon shot noise)
    nor=max(real(S(:)));
    if I==0     %ideal case
        pho=1000000;
        S_loss=zeros(length(S(1,:)));
        for ii=1:length(S(1,:))
            for jj=1:length(S(1,:))
                
                aa=real(S(ii,jj));
                aa=round(pho*aa/nor);
                bb=imag(S(ii,jj));
                bb=round(pho*bb/nor);
                if aa>0
                    Sx=poissrnd(aa);
                elseif aa<0
                    Sx=-poissrnd(-aa);
                end
                if bb>0
                    Sy=poissrnd(bb);
                elseif bb<0
                    Sy=-poissrnd(-bb);
                end
                S_loss(ii,jj)=S(ii,jj)/nor*pho+Sx+1i*Sy;
            end
        end
        S_loss=S_loss./pho.*nor;
    elseif I==1   %photon shot noise with limited amount of photons
        pho=photon;
        S_loss=zeros(length(S(1,:)));
        for ii=1:length(S(1,:))
            for jj=1:length(S(1,:))
                
                aa=real(S(ii,jj));
                aa=round(pho*aa/nor);
                bb=imag(S(ii,jj));
                bb=round(pho*bb/nor);
                if aa>0
                    Sx=poissrnd(aa);
                elseif aa<0
                    Sx=-poissrnd(-aa);
                end
                if bb>0
                    Sy=poissrnd(bb);
                elseif bb<0
                    Sy=-poissrnd(-bb);
                end
                S_loss(ii,jj)=S(ii,jj)/nor*pho+Sx+1i*Sy;
            end
        end
        S_loss=S_loss./pho.*nor;
    end
end