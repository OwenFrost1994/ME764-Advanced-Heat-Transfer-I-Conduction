clear
clc
close all

%% (a)
figure(1)
pdegplot('HW14geometry')
set(gcf,'position',[200,200,1100,1000])
xlabel({'$x[m]$'},'FontSize',20,'Interpreter','Latex');
ylabel({'$y[m]$'},'FontSize',20,'Interpreter','Latex');
set(gca, 'FontName','Times New Roman','FontSize', 20)

%% (b)
%Create a figure showing the mesh using the pdemesh command
[p,e,t]=initmesh('HW14geometry');
[p,e,t]=refinemesh('HW14geometry',p,e,t);
[p,e,t]=refinemesh('HW14geometry',p,e,t);
figure(2)
pdemesh(p,e,t)
set(gcf,'position',[200,200,1100,1000])
xlabel({'$x[m]$'},'FontSize',20,'Interpreter','Latex');
ylabel({'$y[m]$'},'FontSize',20,'Interpreter','Latex');
set(gca, 'FontName','Times New Roman','FontSize', 20)

Rh = 0.02; %m
Rt = 0.005;
W = 0.01;
th = 2.003;
L = 0.025;

kg = 2.5;
kins = 0.25;

h = 100;
tH = 50;
tC = 30;

%point matrix p, edge matrix e and triangle matrix t are given
%number of nodes
nN = size(p,2);
%number of elements
nE = size(t,2);
%number of boundary edges
nB =size(e,2);
en=zeros(3,nB);
en(1,:)=e(1,:);
en(2,:)=e(2,:);
en(3,:)=e(5,:);
e = en;
clear en;

%sd and bc
sd = [kg 0 0 0; kins 0 0 0;]';
bc = [0 0 0; 0 0 0; h tH 0; h tH 0; h tC 0; h tC 0;]';

%Construct the global capacitance and conduction matrix and the generation vector
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
    if bnd<= 6
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
pdeplot(p,e,t,'xydata',tss)
colormap('jet')
colorbar
set(gcf,'position',[200,200,1100,1000])
xlabel({'$x[m]$'},'FontSize',20,'Interpreter','Latex');
ylabel({'$y[m]$'},'FontSize',20,'Interpreter','Latex');
set(gca, 'FontName','Times New Roman','FontSize', 20)
title({'$t_{ss}$ in theground heat exchanger'},'FontSize',20,'Interpreter','Latex')
set(gca, 'FontName','Times New Roman','FontSize', 20)

%% (c)
%the total rate of heat transfer from the hot fluid to the cold
%fluid is the integration of temperature gradient times heat transfer coefficient 
%on the left hole or righthole£º
%q_total = int(kg*dt*ds) on boundary H or boundary C
%or the total heat tranfer on the boundary of hold
%q_total = int(h*(t-tH)*ds) on boundary H or q_total = int(h*(t-tC)*ds) boundary C

[tx,ty] = pdegrad(p,t,tss);%%gradient of t along x and y in center of every element
txn = pdeprtni(p,t,tx);
tyn = pdeprtni(p,t,ty);

figure(4)
quiver(p(1,:),p(2,:),-txn',-tyn')
colormap('jet')
colorbar
set(gcf,'position',[200,200,1100,1000])
xlabel({'$x[m]$'},'FontSize',20,'Interpreter','Latex');
ylabel({'$y[m]$'},'FontSize',20,'Interpreter','Latex');
set(gca, 'FontName','Times New Roman','FontSize', 20)
title({'$t_{ss}$ in theground heat exchanger'},'FontSize',20,'Interpreter','Latex')
set(gca, 'FontName','Times New Roman','FontSize', 20)

qH = 0;
qH1 = 0;
for bn=1:1:nB
    i = e(1,bn);j = e(2,bn);bnd = e(3,bn);
    if bnd == 3 || bnd == 4
        xi = p(1,i); yi = p(2,i);
        xj = p(1,j); yj = p(2,j);
        xij = xj-xi; yij = yj-yi;
        sij = sqrt(xij^2+yij^2);
        he = bc(1,bnd); te = bc(2,bnd); qse = bc(3,bnd);
        qH = qH + sij*kg*(sqrt(txn(i,1)^2+tyn(i,1)^2)+sqrt(txn(j,1)^2+tyn(j,1)^2))/2;
        qH1 = qH1 + sij*he*((tss(i)+tss(j))/2-tH);
    end
end
qH
qH1

qC = 0;
qC1 = 0;
for bn=1:1:nB
    i = e(1,bn);j = e(2,bn);bnd = e(3,bn);
    if bnd == 5 || bnd == 6
        xi = p(1,i); yi = p(2,i);
        xj = p(1,j); yj = p(2,j);
        xij = xj-xi; yij = yj-yi;
        sij = sqrt(xij^2+yij^2);
        he = bc(1,bnd); te = bc(2,bnd); qse = bc(3,bnd);
        qC = qC + sij*kg*(sqrt(txn(i,1)^2+tyn(i,1)^2)+sqrt(txn(j,1)^2+tyn(j,1)^2))/2;
        qC1 = qC1 + sij*he*((tss(i)+tss(j))/2-tC);
    end
end
qC
qC1

% the slight difference between the values comes from the numerical error of
% two different computation method

%% d)
kinsn = [0.1,0.2,0.25,0.5,1.5,2.5,5,10];
thn = 0.001:0.001:0.009;
for kinss=1:1:size(kinsn,2)
    kins = kinsn(1,kinss);
    sd = [kg 0 0 0; kins 0 0 0;]';
    for ths=1:1:size(thn,2)
        [p,e,t]=initmesh(strcat('HW14geometry',num2str(ths)));
        [p,e,t]=refinemesh(strcat('HW14geometry',num2str(ths)),p,e,t);
        [p,e,t]=refinemesh(strcat('HW14geometry',num2str(ths)),p,e,t);
        nN = size(p,2);
        nE = size(t,2);
        nB =size(e,2);
        en=zeros(3,nB);
        en(1,:)=e(1,:);
        en(2,:)=e(2,:);
        en(3,:)=e(5,:);
        e = en;
        clear en;
        
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
                
        qs = zeros(nN,1);
        Ht = zeros(nN,1);
        H = zeros(nN,nN);
        for bn=1:1:nB
            i = e(1,bn);j = e(2,bn);bnd = e(3,bn);
            if bnd<= 6
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
        
        S = K+H;
        r = qs + Ht + g;
        tss = inv(S)*r;
        
        [tx,ty] = pdegrad(p,t,tss);
        txn = pdeprtni(p,t,tx);
        tyn = pdeprtni(p,t,ty);
        qins = 0;
        for bn=1:1:nB
            i = e(1,bn);j = e(2,bn);bnd = e(3,bn);
            if bnd == 7
                xi = p(1,i); yi = p(2,i);
                xj = p(1,j); yj = p(2,j);
                xij = xj-xi; yij = yj-yi;
                sij = sqrt(xij^2+yij^2);
                qins = qins + sij*kins*(sqrt(txn(i,1)^2+tyn(i,1)^2)+sqrt(txn(j,1)^2+tyn(j,1)^2))/2;
            end
        end
        qinsn(kinss,ths) = qins;
    end
end

figure(5)
label = {};
for i = 1:1:size(kinsn,2)
    label = [label,strcat('$k_{ins}$=',num2str(kinsn(i)),'[W/(m*k)]')];
    plot(thn,qinsn(i,:),'LineWidth',2)
    hold on;
end
grid on;
xlim([min(thn) max(thn)])
set(gcf,'position',[200,200,1100,1000])
xlabel({'$th[m]$'},'FontSize',20,'Interpreter','Latex');
ylabel({'$q_{ins}[W/m]$'},'FontSize',20,'Interpreter','Latex');
set(gca, 'FontName','Times New Roman','FontSize', 20)
legend(label,'Interpreter','Latex')

%%It can be seen from the plot that the increase of kins will increase the
%%heat transfer. The reason is that the increase of kins will reduce the
%%resistance on heat transfer through the insulator. And also we can find:
%%When kins < 2.5, the increase of th will reduce the qins. In this case,
%%kg>kins, the insulator works as a barrier of thermal conduction.
%%When kins > 2.5, the qins is a straight line with positive slope. Because
%%kg<kins, increase of th is the increase of total k of heat exchanger.
%%The heat conduction through insulator will be stronger and result in the
%%increase of qins.