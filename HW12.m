clear
clc
close all

dW = 0.5; %m
dH = 0.25;

h = 150;
tinf = 300;

k = 54;
rho = 7600;
c = 810;

gddd = 2000;
qdd = 1000;

%number of nodes
nN = 9;
%number of element
nE = 8;
%number of boundary
nB = 8;

%% a)
%information about the subdomain
sd = [k gddd rho c]';
bc = [h 0 h 0; tinf 0 tinf 0;0 qdd 0 0;];

sd
bc

%% b)
%point matrix contains the coordinates of points
p = [0 0; dW 0; 2*dW 0; dW dH; 0 2*dH; 2*dW 2*dH; dW 3*dH; 0 4*dH; 2*dW 4*dH]';

p

%% c)
%edge matrix contains the nodes of edges and the boundary they belongs to
e = [1 2 1; 2 3 1; 3 6 2; 6 9 2; 9 7 3; 7 8 3; 8 5 4; 5 1 4;]';

e

%% d)
%triangle matrix contains the nodes and the domain
t = [5 1 4 1; 4 1 2 1; 4 2 3 1; 4 3 6 1; 7 4 6 1; 7 5 4 1; 8 5 7 1; 9 7 6 1]';

t

%% e)
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

C
K
g

%% f)
%Construct the specified vectors qs and Ht and the matrix H
qs = zeros(nN,1);
Ht = zeros(nN,1);
H = zeros(nN,nN);
for bn=1:1:nB
    i = e(1,bn);j = e(2,bn);bnd = e(3,bn);
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

qs
Ht
H

%% g)
%Create a figure showing the mesh using the pdemesh command
figure(1)
pdemesh(p,e,t)
xlim([0 2*dW])
ylim([0 4*dH])
set(gcf,'position',[200,300,1000,800])
xlabel({'$x[m]$'},'FontSize',20,'Interpreter','Latex');
ylabel({'$y[m]$'},'FontSize',20,'Interpreter','Latex');
set(gca, 'FontName','Times New Roman','FontSize', 20)

%% g)
%Determine the steady state temperature distribution and create a plot using the pdeplot command
S = K+H;
r = qs + Ht;
tss = S\r;

figure(2)
pdeplot(p,e,t,'xydata',tss)
xlim([0 2*dW])
ylim([0 4*dH])
colorbar
set(gcf,'position',[200,300,1000,800])
xlabel({'$x[m]$'},'FontSize',20,'Interpreter','Latex');
ylabel({'$y[m]$'},'FontSize',20,'Interpreter','Latex');
set(gca, 'FontName','Times New Roman','FontSize', 20)