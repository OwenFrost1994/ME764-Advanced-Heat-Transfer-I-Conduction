clear
clc
close all

D = 0.06; %m
th = 0.008;
L = 0.30;

hinv = 20;
hinl = 300;
tin = 80;

hout = 25;
tout = 300;

k = 10;
rho = 8500;
c = 420;

Davg = D+th;

%% a)
M = 3;
N = 3;
dx = L/(M-1);
dy = pi*Davg/(N-1);

C = Cmake(M,N,rho,c,dx,dy,th);
K = Kmake(M,N,k,dx,dy,th);
Hout = Houtmake(M,N,hout,dx,dy);
Hin = Hinmake(M,N,hinl,hinv,dx,dy);
[tinfin, tinfout] = tinfmake(M,N,tin,tout);
C
K
Hout
Hin
tinfin
tinfout

%% b)
M = 51;
N = 51;
dx = L/(M-1);
dy = pi*Davg/(N-1);

C = Cmake(M,N,rho,c,dx,dy,th);
K = Kmake(M,N,k,dx,dy,th);
Hout = Houtmake(M,N,hout,dx,dy);
Hin = Hinmake(M,N,hinl,hinv,dx,dy);
[tinfin, tinfout] = tinfmake(M,N,tin,tout);

S = K+Hin+Hout;
r = Hin*tinfin+Hout*tinfout;
j=1;
for j=1:1:N
%dirichlet at other x=0 nodes
    dirich(j,1)=M*(j-1)+1;
    dirich(j,2)=tout;
    j=j+1;
end

for n=1:1:size(dirich,1)%loop over every Dirichlet boundary condition
    rw = dirich(n,1);
    for m=1:1:M*N%loop to set all columns of row as 0
        S(rw,m) = 0;
    end
    S(rw,rw) = 1.;%set S(rw,rw) as 1
    r(rw) = dirich(n,2);%set r(rw) as the constant value 
end

tss = S\r;
Tss = reshape(tss,[M,N]);

figure(1)
surfc(0:dx:L,0:dy:pi*Davg,Tss','LineWidth',1)
xlim([0 L])
ylim([0 pi*Davg])
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
%Crank-Nicolson method

dtheta = thetac/2;

%initial condition is tout
tini = tout*ones(M*N,1);

for n=1:1:size(dirich,1)%loop over every Dirichlet boundary condition
    rw = dirich(n,1);
    for m=1:1:M*N%loop to set all columns of row as 0
        C(rw,m) = 0;
    end
    C(rw,rw) = 1.;%set S(rw,rw) as 1
end

A = C+S*dtheta/2;
B = C-S*dtheta/2;
b = r*dtheta;

nstep = 2000;
%transient temperature
tt = zeros(M*N,nstep);
tt(:,1) = tini;

for thetas=1:1:nstep-1
    tt(:,thetas+1) = A\(B*tt(:,thetas)+b);
end

figure(2)
plot(0:dtheta:dtheta*(nstep-1),tt(M,:),'b-','LineWidth',2);
hold on;
plot(0:dtheta:dtheta*(nstep-1),80*ones(1,nstep),'r--','LineWidth',2);
hold on;
plot(0:dtheta:dtheta*(nstep-1),100*ones(1,nstep),'r.-','LineWidth',2);
hold on;
grid on;
xlim([0 dtheta*(nstep-1)])
xlabel({'${\theta}[s]$'},'FontSize',20,'Interpreter','Latex');
ylabel({'$t(x=L,y=0)[K]$'},'FontSize',20,'Interpreter','Latex');
legend({'$t(x=L,y=0)$','$t_{oxegen}$','$t_{ss}$'},'FontSize',20,'Interpreter','Latex');
set(gca, 'FontName','Times New Roman','FontSize', 20)