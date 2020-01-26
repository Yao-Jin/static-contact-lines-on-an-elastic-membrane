close all;
clear;
clc;
% 10/10/2019
% parameterization: X(x,t)=(x,y(x,t)) i.e. a 2-D curve discribed by the height function y(x,t)

% set physical parameters 
cb = 1e-2;
gamma1 = 1.5; % surface tensions on three different surfaces
gamma2 = 1;
gamma3 = 1;
thetas = acos((gamma2-gamma1)/gamma3);
T = 20; % end time
dt = 1e-3;  % modified according to h 
n = 500;  
tol = 1e-8;

% set geometric parameters
xsta = -1.0; xend = 1.0; xlen = 2.0;
yd = 0; yu = 1.0; ylen = 1.0;
kappa0 = 0; kappan = 0;
arearect = 2;
h = xlen/n;
x = (xsta:h:xend)'; % uniformly 

% set initial conditions: a straight line
y = zeros(n+3,1);
% y = zeros(n+1,1);
kappa = zeros(n-1,1);
Radius = 0.4;
xcl1 = -Radius;    ycl1 = yd;
xcl2 = Radius;     ycl2 = yd;
areadrop = Radius^2*pi/2;

pos1 = floor((xcl1-xsta)/h)+1; 
pos2 = floor((xcl2-xsta)/h)+1; % indices impling which line segments our contact pts lie on 
 
itr = 0;    maxitr = T/dt;
yold = y;

% Ypath = zeros(5000,n-1);

while itr<=maxitr
    itr = itr+1; disp(itr);
%     Ypath(itr,:) = y;
    %% Step one: 
    % once know the position of contact pts on the curve
    % using the volume constraint of the droplet to get the thetad
    % the droplet is circle-shaped
    % Xcl1 = (xcl1,ycl1), Xcl2 =(xcl2,ycl2) are on the curve
    
    [thetad,Radius,theta,circenterx,circentery]=reshape_droplet(pos1,xcl1,ycl1,pos2,xcl2,ycl2,x,y,areadrop);
    
    dtheta = abs(thetad-thetas)
    
%     plot(x,y(2:end-1),'-k');
    plot([-1-h;x;1+h],y,'-k');
%     plot(x,y,'-k');
    hold on;
    axis equal;
    para = linspace(theta-pi/2,3/2*pi-theta,200)';
    plot(ones(200,1)*circenterx+Radius*cos(para), ones(200,1)*circentery+Radius*sin(para),'-r');
    plot([xcl1;xcl2],[ycl1;ycl2],'*b','MarkerSize',10);
    legend('membrane','droplet surface','contact points');
    hold off;
    pause(0.001);
    
    %% Step two: 
    % Know previous (kappa, y) to update new (kappa, y)
    % kappa: curvature defined on the nodal pts
    % (x, y(x,t))  position vector 
    [y,kappa] = update_membrane(thetad,Radius,h,xcl1,xcl2,x,y,yu,...
                                kappa,gamma1,gamma2,gamma3,arearect,dt,cb);
                            
    yinfnorm = max(abs(y-yold))
    if yinfnorm<tol 
        save('data_cb1e_2n4000dt1e_3R0_4gamma1_0_5.mat');
        disp('Finished'); 
        break; 
    end
    
    yold = y;
    
    % check if y is symmetric
%     max(abs(y(1:1:end)-y(end:-1:1)))    
    
%     if mod(itr-1,100)==0
%         save('data_cb1e_2n4000dt1e_3R0_4gamma1_0_5.mat');
%     end
    
    %% Step three:
    % Using thetad and new y to update the position of contact pts (xcl(t), y(xcl(t),t))
    [pos1,xcl1,ycl1,pos2,xcl2,ycl2] = move_contactline(pos1,pos2,xcl1,xcl2,gamma3,thetad,thetas,x,y,h,dt,xsta);
    
end
