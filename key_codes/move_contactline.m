function [pos1,xcl1,ycl1,pos2,xcl2,ycl2] = move_contactline(pos1,pos2,xcl1,xcl2,gamma3,thetad,thetas,x,yall,h,dt,xsta)
% Using thetad and new y to update the position of contact pts (xcl(t), y(xcl(t),t))
% idea one: ycl = y(xcl(t),t) update to new xcl then find the corresponding ycl
    y = yall(2:end-1);
    
    Qcl1 = sqrt(1 + ((y(pos1+1)-y(pos1))/h)^2);
    Qcl2 = sqrt(1 + ((y(pos2+1)-y(pos2))/h)^2);
    veloval = gamma3*(cos(thetad)-cos(thetas));
    dist = veloval*dt;
    
    % update contact pts
    xcl1 = xcl1 + dist/Qcl1;
    xcl2 = xcl2 - dist/Qcl2;
    
    pos1 = floor((xcl1-xsta)/h)+1;
    pos2 = floor((xcl2-xsta)/h)+1;
    
    ycl1 = y(pos1) + (y(pos1+1)-y(pos1))*(xcl1-x(pos1))/h; 
    ycl2 = y(pos2) + (y(pos2+1)-y(pos2))*(xcl2-x(pos2))/h;

end

