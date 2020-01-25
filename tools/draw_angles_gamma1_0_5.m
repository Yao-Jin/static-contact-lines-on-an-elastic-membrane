clear;
close all;
clc;
addpath(genpath('.'));
dbstop if error

% gamma1 =1.5; gamma2 = 1; gamma3 = 1;
% C = 0.39564392373896; %>0
% C1 = 2.049145600069344; %>0
% C2 = -2.303604857969651;%<0

gamma1 = 0.5; gamma2 = 1; gamma3 = 1;
C =  0.5352331346596348; %>0
C1 = 2.300744081094803; %>0
C2 = -1.9929416195128653;%<0

phi1 = atan(sqrt(4*gamma1/C^2-1));
phi2 = atan(-sqrt(4*gamma2/C^2-1));
theta0 = pi-2*phi2;

asypo_theta1 = @(s)(theta0-2*(atan(sinh(sqrt(gamma1)*(s+C1)))-phi1)); %s>0
asypo_theta2 = @(s)(theta0-2*(atan(sinh(sqrt(gamma2)*(s+C2)))-phi2)); %s<0


% load("cb1e_4n4000dt1e_3R0_4gamma1_1_5.mat");
load("cb1e_4n4000dt1e_3R0_4gamma1_0_5.mat");
L = 2;
nu = cb/(L^2*gamma3);
xi = sqrt(nu);

hfn = floor(n/2);
linesec = sqrt((y(2:hfn)-y(1:hfn-1)).^2+(x(2:hfn)-x(1:hfn-1)).^2);
s = zeros(hfn,1);
s(1) = 0;
for i = 2:hfn
    s(i) = s(i-1)+linesec(i-1);
end
scl = s(pos1)+sqrt((ycl1-y(pos1))^2+(xcl1-x(pos1))^2);
sbar = (s-scl)/(L*xi);
sbarl = sbar(1:pos1);
sbarr = sbar(pos1+1:hfn);
theta1 = mod(asypo_theta1(sbarr),2*pi);
theta1(theta1(:)<pi) = theta1(theta1(:)<pi)+2*pi;
theta2 = mod(asypo_theta2(sbarl),2*pi);
theta2(theta2(:)<pi) = theta2(theta2(:)<pi)+2*pi;

plot([sbarl;sbarr],[theta2;theta1]-pi,'-b');
hold on;
set(gca,'FontSize',20);
xlabel('$s/\sqrt{\nu}$','Interpreter','latex','FontSize',20);
ylabel('angle','FontSize',20);
% plot([-80 40],[2*pi 2*pi],'k-');
xlim([-10,10]); %-------------------------------------------------------------------


dery = zeros(hfn,1);
dery(1) = (y(2)-y(1))/h;
dery(2:hfn) = (y(3:hfn+1)-y(1:hfn-1))/(2*h);
theta_num = mod(atan(dery),2*pi);
theta_num(theta_num(:)<pi) = theta_num(theta_num(:)<pi)+2*pi;
theta1_num = theta_num(pos1+1:end);
theta2_num = theta_num(1:pos1);
plot(sbar,[theta2_num;theta1_num]-pi,'o-k');


% load("cb1_4e_4n4000dt1e_3R0_4gamma1_1_5.mat");
load("cb1_4e_4n4000dt1e_3R0_4gamma1_0_5.mat");
L = 2;
nu = cb/(L^2*gamma3);
xi = sqrt(nu);

hfn = floor(n/2);
linesec = sqrt((y(2:hfn)-y(1:hfn-1)).^2+(x(2:hfn)-x(1:hfn-1)).^2);
s = zeros(hfn,1);
s(1) = 0;
for i = 2:hfn
    s(i) = s(i-1)+linesec(i-1);
end
scl = s(pos1)+sqrt((ycl1-y(pos1))^2+(xcl1-x(pos1))^2);
sbar = (s-scl)/(L*xi);
dery = zeros(hfn,1);
dery(1) = (y(2)-y(1))/h;
dery(2:hfn) = (y(3:hfn+1)-y(1:hfn-1))/(2*h);
theta_num = mod(atan(dery),2*pi);
theta_num(theta_num(:)<pi) = theta_num(theta_num(:)<pi)+2*pi;
theta1_num = theta_num(pos1+1:end);
theta2_num = theta_num(1:pos1);
plot(sbar,[theta2_num;theta1_num]-pi,'*-k');


% load("cb1_16e_4n4000dt1e_3R0_4gamma1_1_5.mat");
load("cb1_16e_4n4000dt1e_3R0_4gamma1_0_5.mat");
L = 2;
nu = cb/(L^2*gamma3);
xi = sqrt(nu);

hfn = floor(n/2);
linesec = sqrt((y(2:hfn)-y(1:hfn-1)).^2+(x(2:hfn)-x(1:hfn-1)).^2);
s = zeros(hfn,1);
s(1) = 0;
for i = 2:hfn
    s(i) = s(i-1)+linesec(i-1);
end
scl = s(pos1)+sqrt((ycl1-y(pos1))^2+(xcl1-x(pos1))^2);
sbar = (s-scl)/(L*xi);
dery = zeros(hfn,1);
dery(1) = (y(2)-y(1))/h;
dery(2:hfn) = (y(3:hfn+1)-y(1:hfn-1))/(2*h);
theta_num = mod(atan(dery),2*pi);
theta_num(theta_num(:)<pi) = theta_num(theta_num(:)<pi)+2*pi;
theta1_num = theta_num(pos1+1:end);
theta2_num = theta_num(1:pos1);
plot(sbar,[theta2_num;theta1_num]-pi,'s-k');


% load("cb1_64e_4n6000dt1e_3R0_4gamma1_1_5.mat");
load("cb1_64e_4n6000dt1e_3R0_4gamma1_0_5.mat");
L = 2;
nu = cb/(L^2*gamma3);
xi = sqrt(nu);

hfn = floor(n/2);
linesec = sqrt((y(2:hfn)-y(1:hfn-1)).^2+(x(2:hfn)-x(1:hfn-1)).^2);
s = zeros(hfn,1);
s(1) = 0;
for i = 2:hfn
    s(i) = s(i-1)+linesec(i-1);
end
scl = s(pos1)+sqrt((ycl1-y(pos1))^2+(xcl1-x(pos1))^2);
sbar = (s-scl)/(L*xi);
dery = zeros(hfn,1);
dery(1) = (y(2)-y(1))/h;
dery(2:hfn) = (y(3:hfn+1)-y(1:hfn-1))/(2*h);
theta_num = mod(atan(dery),2*pi);
theta_num(theta_num(:)<pi) = theta_num(theta_num(:)<pi)+2*pi;
theta1_num = theta_num(pos1+1:end);
theta2_num = theta_num(1:pos1);
plot(sbar,[theta2_num;theta1_num]-pi,'x-k');

leg = legend({'asymptotic prediction','$c_b = 1\times 10^{-4}$'...
              '$c_b = 2.5\times 10^{-5}$','$c_b = 6.25\times 10^{-6}$'...
              '$c_b = 1.5625\times 10^{-6}$'},'Location','NorthEast');
set(leg,'Interpreter','latex');
set(leg,'FontSize',20);
hold off;

%%---------------------------------------------------------------------------------------------------
axes('position',[0.18,0.19,0.33,0.33]);

% load("cb1e_4n4000dt1e_3R0_4gamma1_1_5.mat");
load("cb1e_4n4000dt1e_3R0_4gamma1_0_5.mat");
hfn = floor(n/2);
linesec = sqrt((y(2:hfn)-y(1:hfn-1)).^2+(x(2:hfn)-x(1:hfn-1)).^2);
s = zeros(hfn,1);
s(1) = 0;
for i = 2:hfn
    s(i) = s(i-1)+linesec(i-1);
end
scl = s(pos1)+sqrt((ycl1-y(pos1))^2+(xcl1-x(pos1))^2);

dery = zeros(hfn,1);
dery(1) = (y(2)-y(1))/h;
dery(2:hfn) = (y(3:hfn+1)-y(1:hfn-1))/(2*h);
theta_num = mod(atan(dery),2*pi);
theta_num(theta_num(:)<pi) = theta_num(theta_num(:)<pi)+2*pi;
theta1_num = theta_num(pos1+1:end);
theta2_num = theta_num(1:pos1);
plot(s-scl,[theta2_num;theta1_num]-pi,'o-k');
hold on;
set(gca,'FontSize',11);
xlabel('$s$','Interpreter','latex','FontSize',20);
ylabel('angle','FontSize',20);
% plot([-80 40],[2*pi 2*pi],'k-');
xlim([-0.05,0.05]); %-------------------------------------------------------------------
% ylim([4.8-pi,6.4-pi]);

% load("cb1_4e_4n4000dt1e_3R0_4gamma1_1_5.mat");
load("cb1_4e_4n4000dt1e_3R0_4gamma1_0_5.mat");
hfn = floor(n/2);
linesec = sqrt((y(2:hfn)-y(1:hfn-1)).^2+(x(2:hfn)-x(1:hfn-1)).^2);
s = zeros(hfn,1);
s(1) = 0;
for i = 2:hfn
    s(i) = s(i-1)+linesec(i-1);
end
scl = s(pos1)+sqrt((ycl1-y(pos1))^2+(xcl1-x(pos1))^2);
dery = zeros(hfn,1);
dery(1) = (y(2)-y(1))/h;
dery(2:hfn) = (y(3:hfn+1)-y(1:hfn-1))/(2*h);
theta_num = mod(atan(dery),2*pi);
theta_num(theta_num(:)<pi) = theta_num(theta_num(:)<pi)+2*pi;
theta1_num = theta_num(pos1+1:end);
theta2_num = theta_num(1:pos1);
plot(s-scl,[theta2_num;theta1_num]-pi,'*-k');


% load("cb1_16e_4n4000dt1e_3R0_4gamma1_1_5.mat");
load("cb1_16e_4n4000dt1e_3R0_4gamma1_0_5.mat");
hfn = floor(n/2);
linesec = sqrt((y(2:hfn)-y(1:hfn-1)).^2+(x(2:hfn)-x(1:hfn-1)).^2);
s = zeros(hfn,1);
s(1) = 0;
for i = 2:hfn
    s(i) = s(i-1)+linesec(i-1);
end
scl = s(pos1)+sqrt((ycl1-y(pos1))^2+(xcl1-x(pos1))^2);
dery = zeros(hfn,1);
dery(1) = (y(2)-y(1))/h;
dery(2:hfn) = (y(3:hfn+1)-y(1:hfn-1))/(2*h);
theta_num = mod(atan(dery),2*pi);
theta_num(theta_num(:)<pi) = theta_num(theta_num(:)<pi)+2*pi;
theta1_num = theta_num(pos1+1:end);
theta2_num = theta_num(1:pos1);
plot(s-scl,[theta2_num;theta1_num]-pi,'s-k');


% load("cb1_64e_4n6000dt1e_3R0_4gamma1_1_5.mat");
load("cb1_64e_4n6000dt1e_3R0_4gamma1_0_5.mat");
hfn = floor(n/2);
linesec = sqrt((y(2:hfn)-y(1:hfn-1)).^2+(x(2:hfn)-x(1:hfn-1)).^2);
s = zeros(hfn,1);
s(1) = 0;
for i = 2:hfn
    s(i) = s(i-1)+linesec(i-1);
end
scl = s(pos1)+sqrt((ycl1-y(pos1))^2+(xcl1-x(pos1))^2);
dery = zeros(hfn,1);
dery(1) = (y(2)-y(1))/h;
dery(2:hfn) = (y(3:hfn+1)-y(1:hfn-1))/(2*h);
theta_num = mod(atan(dery),2*pi);
theta_num(theta_num(:)<pi) = theta_num(theta_num(:)<pi)+2*pi;
theta1_num = theta_num(pos1+1:end);
theta2_num = theta_num(1:pos1);
plot(s-scl,[theta2_num;theta1_num]-pi,'x-k');

hold off;