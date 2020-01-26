function [y,kappa] = update_membrane(thetad,Radius,h,xcl1,xcl2,x,y,yu,kappa,gamma1,gamma2,gamma3,arearect,dt,cb)


    n = length(y)-3;

   % Know previous (kappa, y) to update new (kappa, y)
    % kappa: curvature defined on the nodal pts
    % (x, y(x,t))  position vector 
    % x: with index 1,...,n+1
    % y: with index 1,2,...,n+2,n+3
    % kappa with index 1,2,...,n+2,n+3
    
    % x(1,...,n+1) ---> y(2,...,n+2)
    % x(1,...,n+1) ---> kappa(2,...,n+2) with kappa(2) = kappa(n+2) = 0
    % y(1) and y(n+3) are artificial points
    
    % (Nvecx,Nvecy): outward normal vectors
    % dery, derkappa: spatial derivatives of y & kappa respect to x 
    % der2y, der2kappa: second order spatial derivatives of y & kappa respect to x
    % Q = |d(x,y(x,t))/dx| local length 
    % approximate deforming curves with piecewise linear curves
    dery = zeros(n+1,1);
    dery(1:n+1) = (y(3:n+3)-y(1:n+1))/(2*h);
    
    Q = sqrt(1+dery.^2);
    
    der2y = zeros(n+1,1);
    der2y(1:n+1) = (y(3:n+3)-2*y(2:n+2)+y(1:n+1))/h^2;
    
    % Deal with delta forces by the smoothed delta function
    % ucl1,ucl2 are two parameters characterizing positions of contact line in the domain of parameters
    % Fcl1, Fcl2 are values of delta forces determined by thetad
    % f: force on the whole curve from the delta forces
    Fcl1 = gamma3*sin(thetad);    Fcl2 = gamma3*sin(thetad);
    f = Fcl1./Q.*deltah(x-xcl1*ones(size(x)),h) + Fcl2./Q.*deltah(x-xcl2*ones(size(x)),h);
	
    % piecewise function gamma on the membrane
    gamma = zeros(n+1,1);
    gamma(1:n/2) = (gamma1-gamma2)*steph(x(1:n/2)-xcl1*ones(n/2,1),h)+gamma2*ones(n/2,1);
    gamma(n/2+1:n+1) = (gamma1-gamma2)*steph(xcl2*ones(n/2+1,1)-x(n/2+1:n+1),h)+gamma2*ones(n/2+1,1);
    
    % piecewise function dlammda on the membrane
    % dlemmda = lemmda2-lemmda1 = -kappa3*gamma3 = -1/Radius*gamma3
    dlemmda = zeros(n+1,1);
    dlemmda(1:n/2) = -1/Radius*gamma3*steph(x(1:n/2)-xcl1*ones(n/2,1),h);
    dlemmda(n/2+1:n+1) = -1/Radius*gamma3*steph(xcl2*ones(n/2+1,1)-x(n/2+1:n+1),h);
    
    % update new raw y(x,t) with the semi-implicit scheme
    a = 1./(Q.^2)/h^2+1./(Q.^4).*(dery.*der2y)/(2*h);
    b = -2./(Q.^2)/h^2;
    c = 1./(Q.^2)/h^2-1./(Q.^4).*(dery.*der2y)/(2*h);
    A0 = circshift(diag([0;a(2:end)]),-1)'+diag(b)+circshift(diag([c(1:end-1);0]),1)'; 

    A = dt*diag(Q)*(cb*(A0+diag(0.5*[0;kappa;0].^2))-diag(gamma));
    
    aa = -1./(Q.^3)/h^2;
    bb = 2./(Q.^3)/h^2;
    cc = -1./(Q.^3)/h^2;
    B = circshift(diag([0;aa(2:end)]),-1)'+diag(bb)+circshift(diag([cc(1:end-1);0]),1)'; 
    B = [[aa(1);zeros(n,1)] B [zeros(n,1);cc(end)]];
    
    Stiffmat = [zeros(n+1,1) diag(ones(n+1,1)) zeros(n+1,1) [dt*Q(1)*cb*a(1);zeros(n,1)] A(1:n+1,1:n+1) [zeros(n,1);dt*Q(n+1)*cb*c(end)];
                    B zeros(n+1,1) diag(ones(n+1,1)) zeros(n+1,1);
                    zeros(1,n+4) 1 zeros(1,n+1);
                    zeros(1,2*n+4) 1 0;
                    -gamma2/(2*h) 0 gamma2/(2*h) zeros(1,n) cb/(2*Q(1)*2*h) 0 -cb/(2*Q(1)*2*h) zeros(1,n);
                    zeros(1,n) gamma2/(2*h) 0 -gamma2/(2*h) zeros(1,n) -cb/(2*Q(n+1)*2*h) 0 cb/(2*Q(n+1)*2*h)];
               
    loadvec = zeros(2*(n+3),1);
    loadvec(1:n+1) = y(2:end-1)+dt*(f+dlemmda).*Q;
    
    sol = Stiffmat\loadvec;
    y_raw = sol(2:n+2);
    
%     plot([-1-h;x;1+h],y_rawall,'-o');

    % modification of raw y(x,t) to satisfy the area constraint
    delta = Stiffmat\[dt*Q;zeros(n+5,1)];
    deltay = delta(2:n+2);
    
%     plot(x,deltay,'-o');
    
    areanow = h/2*sum(2*yu*ones(n,1)-y_raw(1:n)-y_raw(2:n+1));
    darea = arearect-areanow;
    lembda2 = darea/(h/2*sum(deltay(1:n)+deltay(2:n+1)));
    y_real = y_raw-lembda2*deltay; 
    
    areanow = h/2*sum(2*yu*ones(n,1)-y_real(1:n)-y_real(2:n+1));
    if abs(areanow-arearect)<1e-10
        disp('Area constraint is satisfied');
    end
    
    % use the modified y to get the new kappa explicitly
    dery = 1/(2*h)*(y_real(3:n+1)-y_real(1:n-1));
    Q = sqrt(1+dery.^2);
    der2y = (y_real(1:n-1)-2*y_real(2:n)+y_real(3:n+1))/h^2;
    kappa = der2y./(Q.^3);
    
    % use kappa(0) = kappa(n) = 0 to get the y(1){y_real(0)} and y(n+3){y_real(n+2)}
    y(2:n+2) = y_real;
    y(1) = 2*y(2)-y(3);
    y(n+3) = 2*y(n+2)-y(n+1);
    
end

%% tool functions:

function result = phih(r)
% separated part for each dimension constructing the whole
% deltah function
result = zeros(size(r));
for j = 1:1:length(r)
    if abs(r(j))<2 
        result(j) = 1/4*(1+cos(pi/2*r(j)));
    end
end
end

function result = deltah(u,h)
% calculate the value of smoothly approximated delta function
% using the form in Peskin(2002)
result = 1/h*phih(u/h);
end

function result = steph(u,h)
% smoothly aoproximated step function
% modified by the integration of deltah
    result = zeros(size(u));
    for j = 1:1:length(u)
        r = u(j)/h;
        if abs(r)<2 
            result(j) = 1/4*(r+2/pi*sin(pi/2*r))+0.5;
        end
        if r>=2
            result(j) = 1;
        end
    end
end