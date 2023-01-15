clear
clc
close all

%% b)
Le = 0.01; %m
Lc = 0.03;
ro = 0.02;
rin = 0.01;

dr = (ro-rin)/2;
dz = (Le+Lc)/2;

P = 10E6;
w = 10;
mu = 0.6;
hb = 150*w;
k = 6.4;
c = 468;
rho = 8238;
tinf = 20+273;
tc = 20+273;
Rc = 0.03;

%%matrix C
C = zeros(9,9);
C(1,1) = rho*c*((rin+dr/2)^2-(rin)^2)*pi*dz/2;
C(2,2) = rho*c*((rin+3*dr/2)^2-(rin+dr/2)^2)*pi*dz/2;
C(3,3) = rho*c*((ro)^2-(rin+3*dr/2)^2)*pi*dz/2;
C(4,4) = 2*C(1,1);
C(5,5) = 2*C(2,2);
C(6,6) = 2*C(3,3);
C(7,7) = C(1,1);
C(8,8) = C(2,2);
C(9,9) = C(3,3);

C
%%matrix K
K = zeros(9,9);
K(1,4) = -k*pi*((rin+dr/2)^2-rin^2)/dz;
K(1,2) = -k*2*pi*(rin+dr/2)*dz/2/dr;
K(1,1) = -(K(1,2)+K(1,4));

K(2,1) = K(1,2);
K(2,5) = -k*pi*((rin+3*dr/2)^2-(rin+dr/2)^2)/dz;
K(2,3) = -k*2*pi*(rin+3*dr/2)*dz/2/dr;
K(2,2) = -(K(2,1)+K(2,5)+K(2,3));

K(3,2) = K(2,3);
K(3,6) = -k*pi*(ro^2-(rin+3*dr/2)^2)/dz;
K(3,3) = -(K(3,6)+K(3,2));

K(4,1) = K(1,4);
K(4,7) = K(1,4);
K(4,5) = 2*K(1,2);
K(4,4) = -(K(4,1)+K(4,5)+K(4,7));

K(5,4) = K(4,5);
K(5,2) = K(2,5);
K(5,8) = K(2,5);
K(5,6) = 2*K(2,3);
K(5,5) = -(K(5,2)+K(5,4)+K(5,6)+K(5,8));

K(6,5) = K(5,6);
K(6,3) = K(3,6);
K(6,9) = -k*pi*(ro^2-(rin+3*dr/2)^2)/dz;
K(6,6) = -(K(6,3)+K(6,5)+K(6,9));

K(7,4) = K(4,7);
K(7,8) = K(1,2);
K(7,7) = -(K(7,4)+K(7,8));

K(8,5) = K(5,8);
K(8,7) = K(7,8);
K(8,9) = K(2,3);
K(8,8) = -(K(8,5)+K(8,7)+K(8,9));

K(9,8) = K(8,9);
K(9,6) = K(6,9);
K(9,9) = -(K(9,6)+K(9,8));

K
%%matrix H
H = zeros(9,9);
H(9,9) = hb*2*pi*ro*dz/2;
H(6,6) = 2*pi*ro*dz/Rc;
H(3,3) = 2*pi*ro*dz/2/Rc;

H
%%vector Q
Q = zeros(9,1);
Q(7) = 0.5*mu*P*2*pi*w*((rin+dr/2)^3-(rin)^3)/3;
Q(8) = 0.5*mu*P*2*pi*w*((rin+3*dr/2)^3-(rin+dr/2)^3)/3;
Q(9) = 0.5*mu*P*2*pi*w*((ro)^3-(rin+3*dr/2)^3)/3;

Q
%%vector Tinf
Tinf = zeros(9,1);
Tinf(3) = tc;
Tinf(6) = tc;
Tinf(9) = tinf;

Tinf


%% c)
S = K+H;
r = H*Tinf+Q;

tss = inv(S)*r;

tss = tss-273;

tss

Tss = reshape(tss,[3,3])';

r = rin:dr:ro;
z = 0:dz:Le+Lc;
[R,Z] = meshgrid(r,z);

rn = rin:dr/20:ro;
zn = 0:dz/100:Le+Lc;
[Rn,Zn] = meshgrid(rn,zn);

Tssn = interp2(R,Z,Tss,Rn,Zn,'spline');

figure(1)
contourf(Rn,Zn,Tssn,'LineWidth',1)
xlim([rin ro])
ylim([0 Le+Lc])
colorbar
set(gcf,'position',[200,300,1000,800])
xlabel({'$r[m]$'},'FontSize',20,'Interpreter','Latex');
ylabel({'$z[m]$'},'FontSize',20,'Interpreter','Latex');
set(gca, 'FontName','Times New Roman','FontSize', 20)

figure(2)
surf(Rn,Zn,Tssn)
xlim([rin ro])
ylim([0 Le+Lc])
colorbar
set(gcf,'position',[200,300,1000,800])
xlabel({'$r[m]$'},'FontSize',20,'Interpreter','Latex');
ylabel({'$z[m]$'},'FontSize',20,'Interpreter','Latex');
zlabel({'$t[{}^{\circ}C]$'},'FontSize',20,'Interpreter','Latex');
set(gca, 'FontName','Times New Roman','FontSize', 20)