function [thetad,Radius,theta,circenterx,circentery]=reshape_droplet(pos1,xcl1,ycl1,pos2,xcl2,ycl2,x,yall,areadrop)

% once know the position of contact pts on the curve
% using the volume constraint of the droplet to get the thetad
% the droplet is circle-shaped
% Xcl1 = (xcl1,ycl1), Xcl2 =(xcl2,ycl2) are on the curve
    y = yall(2:end-1);
    
    if mod(pos1+1+pos2,2)==0 
%         xmid = x((pos1+1+pos2)/2); 
        ymid = y((pos1+1+pos2)/2);
    else
        left = (pos1+pos2)/2; right = left + 1;
%         xmid = (x(left)+x(right))/2;    
        ymid = (y(left)+y(right))/2;
    end
    
    xc = (xcl1+xcl2)/2;
    yc = (ycl1+ycl2)/2; % ycl1 == ycl2 symmetric
    L = sqrt((xcl2-xcl1)^2+(ycl2-ycl1)^2)/2;

    if  yc>=ymid  % the sign above "+" concave;   "-" convex 
        darea = areadrop - polyarea([xcl1;x(pos1+1:pos2);xcl2], [ycl1;y(pos1+1:pos2);ycl2]);
    else
        darea = areadrop + polyarea([xcl1;x(pos1+1:pos2);xcl2], [ycl1;y(pos1+1:pos2);ycl2]);
    end 
    
    
    % Use dichotomy to find theta_d
    lbd = 0.001*pi; rbd = 0.999*pi;  % set a proper initial range of theta for dichotomy
    tol = 1e-9;
    while (rbd-lbd>tol)
        theta = 0.5*(lbd+rbd);
        R = L/sin(theta);
        areanow = (pi-theta)*R^2+L^2/tan(theta);
        if areanow < darea 
            rbd = theta;
        else
            lbd = theta;
        end
    end

   % Once find the proper theta, compute the thetad
    circenterx = xc;    circentery = yc + L/tan(theta);
    R = L/sin(theta);
    Radius = R; % update the new Radius (previous radius)
    x1 = circenterx+R*cos((3/2*pi-theta)-1e-3);
    y1 = circentery+R*sin((3/2*pi-theta)-1e-3);
    clvecx = x1-xcl1;
    clvecy = y1-ycl1;
    cllen = sqrt(clvecx^2+clvecy^2);
    
    vecxl = x(pos1+1)-x(pos1-1);
    vecyl = y(pos1+1)-y(pos1-1);
    
    vecxr = x(pos1+2)-x(pos1);
    vecyr = y(pos1+2)-y(pos1);
    
    vecx = (x(pos1+1)-xcl1)/(x(pos1+1)-x(pos1))*vecxl+(xcl1-x(pos1))/(x(pos1+1)-x(pos1))*vecxr;
    vecy = (x(pos1+1)-xcl1)/(x(pos1+1)-x(pos1))*vecyl+(xcl1-x(pos1))/(x(pos1+1)-x(pos1))*vecyr;
    len = sqrt(vecx^2+vecy^2);
    
    thetad = acos((vecx*clvecx+vecy*clvecy)/(len*cllen));
    
%     % Display the figure of the current droplet
%     plot([xl; x; xr],[ysta; y; yend],'-k');
%     hold on;
%     %plot([0],[0],'.','MarkerSize',13);
%     plot([xcl1;xcl2],[ycl1;ycl2],'*');
%     %axis([xl xr yd yu]);
%     axis equal;
%     para = linspace(theta-pi/2,3/2*pi-theta,100)';
%     plot(ones(100,1)*circenterx+R*cos(para), ones(100,1)*circentery+R*sin(para),'-r');
%     hold off;
%     %pause(0.01);


end