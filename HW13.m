clc
close all

W = 0.015; %m
L = 0.015;

km = 11.6;
cm = 320;
rhom = 6500;
kf = 0.25;
cf = 3100;
rhof = 1100;

qdd = 1000;
h = 1E10;
tinf = 100;

%point matrix p, edge matrix e and triangle matrix t are given by pdetool box
%then modify them
%number of nodes
nN = size(p_pdetool,2);
p=p_pdetool;
%number of elements
nE = size(t_pdetool,2);
t=t_pdetool;
%number of boundary edges
nB =size(e_pdetool,2);
e=zeros(3,nB);
e(1,:)=e_pdetool(1,:);
e(2,:)=e_pdetool(2,:);
e(3,:)=e_pdetool(5,:);

%% effective conductivity in the x-direction
%sd and bc
sd = [kf 0 rhof cf; km 0 rhom cm; km 0 rhom cm]';
bc = [0 0 0; 0 0 0; h tinf 0; h tinf 0; 0 0 0; 0 0 0; 0 0 qdd; 0 0 qdd;]';

%Construct the global capacitance and conduction matrix and the generation vector
C = zeros(nN,nN);
K = zeros(nN,nN);
g = zeros(nN,1);
for en=1:1:nE
    i = t(1,en);j = t(2,en);k = t(3,en);d = t(4,en);
    xi = p(1,i); yi = p(2,i);
    xj = p(1,j); yj = p(2,j);
    xk = p(1,k); yk = p(2,k);
    xij = xj-xi; yij = yj-yi;
    xik = xk-xi; yik = yk-yi;
    xjk = xk-xj; yjk = yk-yj;
    bijk = xij*yjk-xjk*yij;
    Ae = bijk/2;
    ke = sd(1,d); ge = sd(2,d); rhoe = sd(3,d); ce = sd(4,d);
    C(i,i) = C(i,i) + rhoe*ce*Ae/6;
    C(i,j) = C(i,j) + rhoe*ce*Ae/12;
    C(i,k) = C(i,k) + rhoe*ce*Ae/12;
    C(j,i) = C(j,i) + rhoe*ce*Ae/12;
    C(j,j) = C(j,j) + rhoe*ce*Ae/6;
    C(j,k) = C(j,k) + rhoe*ce*Ae/12;
    C(k,i) = C(k,i) + rhoe*ce*Ae/12;
    C(k,j) = C(k,j) + rhoe*ce*Ae/12;
    C(k,k) = C(k,k) + rhoe*ce*Ae/6;
    
    K(i,i) = K(i,i) + ke*Ae*(yjk^2+xjk^2)/bijk^2;
    K(i,j) = K(i,j) + ke*Ae*(-yik*yjk-xik*xjk)/bijk^2;
    K(i,k) = K(i,k) + ke*Ae*(yjk*yij+xjk*xij)/bijk^2;
    K(j,i) = K(j,i) + ke*Ae*(-yik*yjk-xik*xjk)/bijk^2;
    K(j,j) = K(j,j) + ke*Ae*(yik^2+xik^2)/bijk^2;
    K(j,k) = K(j,k) + ke*Ae*(-yik*yij-xik*xij)/bijk^2;
    K(k,i) = K(k,i) + ke*Ae*(yjk*yij+xjk*xij)/bijk^2;
    K(k,j) = K(k,j) + ke*Ae*(-yik*yij-xik*xij)/bijk^2;
    K(k,k) = K(k,k) + ke*Ae*(yij^2+xij^2)/bijk^2;
    
    g(i,1) = g(i,1) + ge*Ae/3;
    g(j,1) = g(j,1) + ge*Ae/3;
    g(k,1) = g(k,1) + ge*Ae/3;
end

%Construct the specified vectors qs and Ht and the matrix H
qs = zeros(nN,1);
Ht = zeros(nN,1);
H = zeros(nN,nN);
for bn=1:1:nB
    i = e(1,bn);j = e(2,bn);bnd = e(3,bn);
    if bnd<= 8
        xi = p(1,i); yi = p(2,i);
        xj = p(1,j); yj = p(2,j);
        xij = xj-xi; yij = yj-yi;
        sij = sqrt(xij^2+yij^2);
        he = bc(1,bnd); te = bc(2,bnd); qse = bc(3,bnd);
        qs(i,1) = qs(i,1) + qse*sij/2;
        qs(j,1) = qs(j,1) + qse*sij/2;
        
        Ht(i,1) = Ht(i,1) + he*te*sij/2;
        Ht(j,1) = Ht(j,1) + he*te*sij/2;
        
        H(i,i) = H(i,i) + he*sij/3;
        H(i,j) = H(i,j) + he*sij/6;
        H(j,i) = H(j,i) + he*sij/6;
        H(j,j) = H(j,j) + he*sij/3;
    end
end

%Create a figure showing the mesh using the pdemesh command
figure(1)
pdemesh(p,e,t)
xlim([0 W])
ylim([0 L])
set(gcf,'position',[200,300,1000,800])
xlabel({'$x[m]$'},'FontSize',20,'Interpreter','Latex');
ylabel({'$y[m]$'},'FontSize',20,'Interpreter','Latex');
set(gca, 'FontName','Times New Roman','FontSize', 20)

%Determine the steady state temperature distribution
S = K+H;
r = qs + Ht + g;
tss = inv(S)*r;

%Create a plot using the pdeplot command
figure(2)
pdesurf(p,t,tss)
colormap('jet')
xlim([0 W])
ylim([0 L])
colorbar
set(gcf,'position',[200,300,1000,800])
xlabel({'$x[m]$'},'FontSize',20,'Interpreter','Latex');
ylabel({'$y[m]$'},'FontSize',20,'Interpreter','Latex');
title({'$t_{ss}$ for effective conductivity in x'},'FontSize',20,'Interpreter','Latex')
set(gca, 'FontName','Times New Roman','FontSize', 20)

tavex = 0;
nnx = 0;
for bn=1:1:nB
    if e(3,bn)==7 || e(3,bn)==8
        tavex = tavex + (tss(e(1,bn))+tss(e(2,bn)))/2;
        nnx = nnx+1;
    end
end
tavex=tavex/nnx;

keffx = qdd*L/(tavex - tinf);

%% effective conductivity in the x-direction is
fprintf('Effective conductivity in the x-direction is k_{effx} %f [W/m-K]\n', keffx);

%% effective conductivity in the y-direction
%sd and bc
sd = [kf 0 rhof cf; km 0 rhom cm; km 0 rhom cm]';
bc = [h tinf 0; h tinf 0; 0 0 0; 0 0 0; 0 0 qdd; 0 0 qdd; 0 0 0; 0 0 0;]';

%Construct the global capacitance and conduction matrix and the generation vector
C = zeros(nN,nN);
K = zeros(nN,nN);
g = zeros(nN,1);
for en=1:1:nE
    i = t(1,en);j = t(2,en);k = t(3,en);d = t(4,en);
    xi = p(1,i); yi = p(2,i);
    xj = p(1,j); yj = p(2,j);
    xk = p(1,k); yk = p(2,k);
    xij = xj-xi; yij = yj-yi;
    xik = xk-xi; yik = yk-yi;
    xjk = xk-xj; yjk = yk-yj;
    bijk = xij*yjk-xjk*yij;
    Ae = bijk/2;
    ke = sd(1,d); ge = sd(2,d); rhoe = sd(3,d); ce = sd(4,d);
    C(i,i) = C(i,i) + rhoe*ce*Ae/6;
    C(i,j) = C(i,j) + rhoe*ce*Ae/12;
    C(i,k) = C(i,k) + rhoe*ce*Ae/12;
    C(j,i) = C(j,i) + rhoe*ce*Ae/12;
    C(j,j) = C(j,j) + rhoe*ce*Ae/6;
    C(j,k) = C(j,k) + rhoe*ce*Ae/12;
    C(k,i) = C(k,i) + rhoe*ce*Ae/12;
    C(k,j) = C(k,j) + rhoe*ce*Ae/12;
    C(k,k) = C(k,k) + rhoe*ce*Ae/6;
    
    K(i,i) = K(i,i) + ke*Ae*(yjk^2+xjk^2)/bijk^2;
    K(i,j) = K(i,j) + ke*Ae*(-yik*yjk-xik*xjk)/bijk^2;
    K(i,k) = K(i,k) + ke*Ae*(yjk*yij+xjk*xij)/bijk^2;
    K(j,i) = K(j,i) + ke*Ae*(-yik*yjk-xik*xjk)/bijk^2;
    K(j,j) = K(j,j) + ke*Ae*(yik^2+xik^2)/bijk^2;
    K(j,k) = K(j,k) + ke*Ae*(-yik*yij-xik*xij)/bijk^2;
    K(k,i) = K(k,i) + ke*Ae*(yjk*yij+xjk*xij)/bijk^2;
    K(k,j) = K(k,j) + ke*Ae*(-yik*yij-xik*xij)/bijk^2;
    K(k,k) = K(k,k) + ke*Ae*(yij^2+xij^2)/bijk^2;
    
    g(i,1) = g(i,1) + ge*Ae/3;
    g(j,1) = g(j,1) + ge*Ae/3;
    g(k,1) = g(k,1) + ge*Ae/3;
end

%Construct the specified vectors qs and Ht and the matrix H
qs = zeros(nN,1);
Ht = zeros(nN,1);
H = zeros(nN,nN);
for bn=1:1:nB
    i = e(1,bn);j = e(2,bn);bnd = e(3,bn);
    if bnd<= 8
        xi = p(1,i); yi = p(2,i);
        xj = p(1,j); yj = p(2,j);
        xij = xj-xi; yij = yj-yi;
        sij = sqrt(xij^2+yij^2);
        he = bc(1,bnd); te = bc(2,bnd); qse = bc(3,bnd);
        qs(i,1) = qs(i,1) + qse*sij/2;
        qs(j,1) = qs(j,1) + qse*sij/2;
        
        Ht(i,1) = Ht(i,1) + he*te*sij/2;
        Ht(j,1) = Ht(j,1) + he*te*sij/2;
        
        H(i,i) = H(i,i) + he*sij/3;
        H(i,j) = H(i,j) + he*sij/6;
        H(j,i) = H(j,i) + he*sij/6;
        H(j,j) = H(j,j) + he*sij/3;
    end
end

%Determine the steady state temperature distribution
S = K+H;
r = qs + Ht + g;
tss = inv(S)*r;

%Create a plot using the pdeplot command
figure(3)
pdesurf(p,t,tss)
colormap('jet')
xlim([0 W])
ylim([0 L])
colorbar
set(gcf,'position',[200,300,1000,800])
xlabel({'$x[m]$'},'FontSize',20,'Interpreter','Latex');
ylabel({'$y[m]$'},'FontSize',20,'Interpreter','Latex');
title({'$t_{ss}$ for effective conductivity in y'},'FontSize',20,'Interpreter','Latex')
set(gca, 'FontName','Times New Roman','FontSize', 20)

tavex = 0;
nnx = 0;
for bn=1:1:nB
    if e(3,bn)==5 || e(3,bn)==6
        tavex = tavex + (tss(e(1,bn))+tss(e(2,bn)))/2;
        nnx = nnx+1;
    end
end
tavex=tavex/nnx;

keffy = qdd*W/(tavex - tinf);

%% effective conductivity in the y-direction is
fprintf('Effective conductivity in the x-direction is k_{effy} %f [W/m-K]\n', keffy);