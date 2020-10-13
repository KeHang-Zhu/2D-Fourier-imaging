function [z,B0]=recon_2D(image,kx,ky)
    px=1000/kx;
    py=1000/ky;
    val=image(:,3);
    num=length(val(:));
    x=image(:,1);
    y=image(:,2);
    
%% move all the points to one small region
    for ii=1:num
        if (image(ii,1)>px)
            while (x>px)
                x=x-px;
            end
        end
        if (image(ii,2)>py)
            while (y>px)
                y=y-py;
            end
        end
    end
    [B0,num0]=max(val(:)); % set the initial value for fitting B0
    x0=x(num0);            % set the initial value for fitting x0
    y0=y(num0);            % set the initial value for fitting y0
    
%% plot the small region 
    figure
    plot3D(x,y,val);
    xlabel('um')
    ylabel('um')
    zlabel('magnetic field strength')
    
%% make some modification
    
    
    
end