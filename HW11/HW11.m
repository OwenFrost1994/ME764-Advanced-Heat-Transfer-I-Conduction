clear
clc
close all

%% Part1
W = 0.005; %m

h = 1100;
tinf = 4.2;

k = 10;
rho = 8000;
c = 120;

gdd = 1E5;

%% a)
M = 3;
N = 3;
dx = W/(M-1);
dy = W/(N-1);

C = Cmake(M,N,rho,c,dx,dy);
K = Kmake(M,N,k,dx,dy);
H = Hmake(M,N,h,dx,dy);
Tinf = tinfmake(M,N,tinf);
g = gmake(M,N,dx,dy,gdd);

C
K
H
Tinf
g

%% b)
M = 21;
N = 21;
dx = W/(M-1);
dy = W/(N-1);

C = Cmake(M,N,rho,c,dx,dy);
K = Kmake(M,N,k,dx,dy);
H = Hmake(M,N,h,dx,dy);
Tinf = tinfmake(M,N,tinf);
g = gmake(M,N,dx,dy,gdd);

S = K+H;
r = H*Tinf+g;

tss = S\r;
Tss = reshape(tss,[M,N]);

figure(1)
surfc(0:dx:W,0:dy:W,Tss','LineWidth',1)
xlim([0 W])
ylim([0 W])
colorbar
set(gcf,'position',[200,300,1000,800])
xlabel({'$x[m]$'},'FontSize',20,'Interpreter','Latex');
ylabel({'$y[m]$'},'FontSize',20,'Interpreter','Latex');
zlabel({'$t_{ss}[K]$'},'FontSize',20,'Interpreter','Latex');
set(gca, 'FontName','Times New Roman','FontSize', 20)

%% c)
%the critical time should be
thetac = 2/max(sum(abs(S),2)./diag(C));
thetac

%% d)
%the equilibrium time should be
thetaeq = W*2.197/(h/(rho*c));
thetaeq

%% e)
%Crank-Nicolson method

dtheta = thetac/2;

%initial condition is
tini = tinf*ones(M*N,1);

A = C+S*dtheta/2;
B = C-S*dtheta/2;
b = r*dtheta;

nstep = floor(15/dtheta)+1;
%transient temperature
tt = zeros(M*N,nstep);
tt(:,1) = tini;

for thetas=1:1:nstep-1
    tt(:,thetas+1) = A\(B*tt(:,thetas)+b);
end

figure(2)
plot(0:dtheta:dtheta*(nstep-1),tt(M,:),'b-','LineWidth',2);
hold on;
plot(0:dtheta:dtheta*(nstep-1),tt(1,:),'r--','LineWidth',2);
hold on;
plot(0:dtheta:dtheta*(nstep-1),tt(M*(N-1)+1,:),'k.-','LineWidth',2);
hold on;
grid on;
set(gcf,'position',[200,300,1000,800])
xlim([0 dtheta*(nstep-1)])
xlabel({'${\theta}[s]$'},'FontSize',20,'Interpreter','Latex');
ylabel({'$t[K]$'},'FontSize',20,'Interpreter','Latex');
legend({'$t(x=W,y=0)$','$t(x=0,y=0)$','$t(x=0,y=W)$'},'Location','southeast','FontSize',20,'Interpreter','Latex');
set(gca, 'FontName','Times New Roman','FontSize', 20)

%% f)
%exact solution
order = 6;

[X,lambda] = eigs(S,C,order,'SM');

for i=1:1:order
    b(i) = X(:,i)'*C*(tini-tss);
end

%transient temperature
te = zeros(M*N,nstep);

for thetas=1:1:nstep
    theta = (thetas-1)*dtheta;
    for j=1:1:order
        te(:,thetas) = te(:,thetas) + b(j)*X(:,j)*exp(-lambda(j,j)*theta);
    end
    te(:,thetas) = te(:,thetas) + tss;
end



figure(3)
plot(0:dtheta:dtheta*(nstep-1),tt(M,:),'b-','LineWidth',2);
hold on;
plot(0:dtheta:dtheta*(nstep-1),tt(1,:),'r-','LineWidth',2);
hold on;
plot(0:dtheta:dtheta*(nstep-1),tt(M*(N-1)+1,:),'k-','LineWidth',2);
hold on;
plot(0:dtheta:dtheta*(nstep-1),te(M,:),'b--','LineWidth',4);
hold on;
plot(0:dtheta:dtheta*(nstep-1),te(1,:),'r--','LineWidth',4);
hold on;
plot(0:dtheta:dtheta*(nstep-1),te(M*(N-1)+1,:),'k--','LineWidth',4);
hold on;
grid on;
set(gcf,'position',[200,300,1000,800])
xlim([0 dtheta*(nstep-1)])
xlabel({'${\theta}[s]$'},'FontSize',20,'Interpreter','Latex');
ylabel({'$t[K]$'},'FontSize',20,'Interpreter','Latex');
legend({'$t(x=W,y=0)$','$t(x=0,y=0)$','$t(x=0,y=W)$','$t_{e}(x=W,y=0)$','$t_{e}(x=0,y=0)$','$t_{e}(x=0,y=W)$'},'Location','southeast','FontSize',20,'Interpreter','Latex');
set(gca, 'FontName','Times New Roman','FontSize', 20)

%% Part2
W = 0.05; %m
L = 0.15;

tinf = 230;
tini = 320;

hfd = 10;
dhdev = 500;
beta = 5;

rho = 3200;

%% a)
M = 31;
N = 31;
dx = W/(M-1);
dy = L/(N-1);

%BC on left side
H = Hymake(M,N,hfd,dhdev,beta,dx,dy,L);
Tinf = tinfset(M,N,tinf);

%Crank-Nicolson method
dtheta = 10;
ttheta = 80000;
nstep = floor(ttheta/dtheta)+1;

%transient temperature
tt = zeros(M*N,nstep);
%initial condition must be defined before the computation
tt(:,1) = tini*ones(M*N,1);

for thetas=1:1:nstep-1
    C = Ctmake(M,N,rho,dx,dy,tt(:,thetas));
    K = Ktmake(M,N,k,dx,dy,tt(:,thetas));
    
    S = K+H;
    r = H*Tinf;
    
    A = C+S*dtheta/2;
    B = C-S*dtheta/2;
    b = r*dtheta;
    
    tt(:,thetas+1) = A\(B*tt(:,thetas)+b);
end

figure(4)
surfc(0:dx:W,0:dy:L,reshape(tt(:,floor(3000/dtheta)+1),[M,N])','LineWidth',1)
xlim([0 W])
ylim([0 L])
colorbar
set(gcf,'position',[200,300,1000,800])
xlabel({'$x[m]$'},'FontSize',20,'Interpreter','Latex');
ylabel({'$y[m]$'},'FontSize',20,'Interpreter','Latex');
zlabel({'$t[K]$'},'FontSize',20,'Interpreter','Latex');
title('t=3000s')
set(gca, 'FontName','Times New Roman','FontSize', 20)
figure(5)
surfc(0:dx:W,0:dy:L,reshape(tt(:,floor(6000/dtheta)+1),[M,N])','LineWidth',1)
xlim([0 W])
ylim([0 L])
colorbar
set(gcf,'position',[200,300,1000,800])
xlabel({'$x[m]$'},'FontSize',20,'Interpreter','Latex');
ylabel({'$y[m]$'},'FontSize',20,'Interpreter','Latex');
zlabel({'$t[K]$'},'FontSize',20,'Interpreter','Latex');
title('t=6000s')
set(gca, 'FontName','Times New Roman','FontSize', 20)
figure(6)
surfc(0:dx:W,0:dy:L,reshape(tt(:,floor(9000/dtheta)+1),[M,N])','LineWidth',1)
xlim([0 W])
ylim([0 L])
colorbar
set(gcf,'position',[200,300,1000,800])
xlabel({'$x[m]$'},'FontSize',20,'Interpreter','Latex');
ylabel({'$y[m]$'},'FontSize',20,'Interpreter','Latex');
zlabel({'$t[K]$'},'FontSize',20,'Interpreter','Latex');
title('t=9000s')
set(gca, 'FontName','Times New Roman','FontSize', 20)

%the temperature after long time
figure(7)
surfc(0:dx:W,0:dy:L,reshape(tt(:,floor(20000/dtheta)+1),[M,N])','LineWidth',1)
xlim([0 W])
ylim([0 L])
colorbar
set(gcf,'position',[200,300,1000,800])
xlabel({'$x[m]$'},'FontSize',20,'Interpreter','Latex');
ylabel({'$y[m]$'},'FontSize',20,'Interpreter','Latex');
zlabel({'$t[K]$'},'FontSize',20,'Interpreter','Latex');
title('t=20000s')
set(gca, 'FontName','Times New Roman','FontSize', 20)
figure(8)
surfc(0:dx:W,0:dy:L,reshape(tt(:,floor(40000/dtheta)+1),[M,N])','LineWidth',1)
xlim([0 W])
ylim([0 L])
colorbar
set(gcf,'position',[200,300,1000,800])
xlabel({'$x[m]$'},'FontSize',20,'Interpreter','Latex');
ylabel({'$y[m]$'},'FontSize',20,'Interpreter','Latex');
zlabel({'$t[K]$'},'FontSize',20,'Interpreter','Latex');
title('t=40000s')
set(gca, 'FontName','Times New Roman','FontSize', 20)
figure(9)
surfc(0:dx:W,0:dy:L,reshape(tt(:,floor(80000/dtheta)+1),[M,N])','LineWidth',1)
xlim([0 W])
ylim([0 L])
colorbar
set(gcf,'position',[200,300,1000,800])
xlabel({'$x[m]$'},'FontSize',20,'Interpreter','Latex');
ylabel({'$y[m]$'},'FontSize',20,'Interpreter','Latex');
zlabel({'$t[K]$'},'FontSize',20,'Interpreter','Latex');
title('t=80000s')
set(gca, 'FontName','Times New Roman','FontSize', 20)
