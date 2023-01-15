clear
clc
close all

W = 0.05;
H = 0.02;
c = 0.005;
thb = 0.0035;
thf = 0.001;
thg = 0.001;
Lf = 0.01;
hb = 32;
Tinf = 20;
qdd = 20E4;
k = 50;
kf = 125;
Rdd = 1E-4;

%% b)
% the intervals of eigenvalue \lambda and the values of \lambda
N = 10;
y = 0:N*pi/(100*N):N*pi;
for i = 1:1:size(y,2)
    res(i) = y(i)*sin(y(i))/W-hb*cos(y(i))/k;
end

%%solve \lambda
fun = @(x) x*sin(x)/W-hb*cos(x)/k;
label = {};
for i = 1:1:N
    label = [label,strcat(num2str(i),'\pi')];
    x0 = (i-1)*pi;
    x1 = 0.5*pi+(i-1)*pi;
    lambda(i) = fzero(fun,[x0,x1])/W;
end
label = ['0',label];

figure(1)
subplot(1,2,1)
plot(y,res,'b-','LineWidth',2);
hold on;
for i = 1:1:N
    plot([0.5*pi+(i-1)*pi,0.5*pi+(i-1)*pi],[min(res),max(res)],'--k','LineWidth',2)
    hold on;
end
plot(lambda*W,zeros(size(lambda)),'.r','Markersize',30,'LineWidth',2);
hold on;
grid on;
xlim([0 N*pi])
set(gca,'xtick',0:pi:10*pi);
set(gca,'xticklabels',label);
xlabel({'$\lambda_{i}th$'},'FontSize',20,'Interpreter','Latex');
ylabel({'$res$'},'FontSize',20,'Interpreter','Latex');
set(gca, 'FontName','Times New Roman','FontSize', 20)
subplot(1,2,2)
plot(y,abs(res),'b-','LineWidth',2);
hold on;
for i = 1:1:N
    plot([0.5*pi+(i-1)*pi,0.5*pi+(i-1)*pi],[min(abs(res)),max(abs(res))],'--k','LineWidth',2)
    hold on;
end
plot(lambda*W,zeros(size(lambda)),'.r','Markersize',30,'LineWidth',2);
hold on;
grid on;
xlim([0 N*pi])
set(gca,'xtick',0:pi:10*pi);
set(gca,'xticklabels',label);
xlabel({'$\lambda_{i}th$'},'FontSize',20,'Interpreter','Latex');
ylabel({'$abs(res)$'},'FontSize',20,'Interpreter','Latex');
set(gca, 'FontName','Times New Roman','FontSize', 20)

%% c)
%%hbeff for computation
Rb = thb/kf;
Rg = (thg+thf)/(hb*thg);
m = sqrt(2*hb/(kf*thf));
Rf = m*(thg+thf)/(2*tanh(m*Lf)*hb);
R = Rdd + Rb + 1/(1/Rf+1/Rg);

hbeff = 1/R;

N = 50;
%%solve \lambda and Ci
fun = @(x) x*sin(x)/W-hb*cos(x)/k;
for i=1:1:N
    x0 = (i-1)*pi;
    x1 = 0.5*pi+(i-1)*pi;
    lambda(i) = fzero(fun,[x0,x1])/W;
    C(i) = qdd*sin(lambda(i)*c)/(lambda(i)*k)/((lambda(i)*sinh(lambda(i)*H)+hbeff*cosh(lambda(i)*H)/k)*(W/2+sin(2*lambda(i)*W)/(4*lambda(i))));
end

Y = 20;
X = 100;
y = 0:H/Y:H;
x = 0:W/X:W;
t = zeros(size(y,2),size(x,2));
for i=1:1:size(y,2)
    for j=1:1:size(x,2)
        for m=1:1:N
            t(i,j) = t(i,j)+cos(lambda(m)*x(j))*(C(m)*cosh(lambda(m)*y(i))+hbeff*C(m)*sinh(lambda(m)*y(i))/(k*lambda(m)));
        end
    end
end
t = t+20;

figure(2)
contourf(0:1/X:1,0:1/Y:1,t,'LineWidth',0)
colorbar
xlabel({'$x/W$'},'FontSize',20,'Interpreter','Latex');
ylabel({'$y/H$'},'FontSize',20,'Interpreter','Latex');
set(gca, 'FontName','Times New Roman','FontSize', 20)
title('Temperature distribution in the spreader plate')

%% d)

%%hbeff for computation
Rb = thb/kf;
Rg = (thg+thf)/(hb*thg);
m = sqrt(2*hb/(kf*thf));
Rf = m*(thg+thf)/(2*tanh(m*Lf)*hb);
R = Rdd + Rb + 1/(1/Rf+1/Rg);

hbeff = 1/R;

N = 50;
Wm = [0.05,0.10,0.2,0.3,0.4,0.5];
Hn = 0.02/100:0.02/100:0.20;
maxt = zeros(size(Wm,2),size(Hn,2));
% for m=1:1:size(Wm,2)
%     W = Wm(m);
%     for n=1:1:size(Hn,2)
%         H = Hn(n);
%         %solve \lambda and Ci
%         fun = @(x) x*sin(x)/W-hb*cos(x)/k;
%         for i=1:1:N
%             x0 = (i-1)*pi;
%             x1 = 0.5*pi+(i-1)*pi;
%             lambda(i) = fzero(fun,[x0,x1])/W;
%             C(i) = qdd*sin(lambda(i)*c)/(lambda(i)*k)/((lambda(i)*sinh(lambda(i)*H)+hbeff*cosh(lambda(i)*H)/k)*(W/2+sin(2*lambda(i)*W)/(4*lambda(i))));
%         end
%         Xn = 0:0.01/20:W;
%         t = zeros(1,size(Xn,2));
%         for i=1:1:size(Xn,2)
%             for j = 1:1:N
%                 t(1,i) = t(1,i)+cos(lambda(j)*Xn(i))*(C(j)*cosh(lambda(j)*H)+hbeff*C(j)*sinh(lambda(j)*H)/(k*lambda(j)));
%             end
%         end
%         maxt(m,n) = max(t);
%     end
% end
for m=1:1:size(Wm,2)
    W = Wm(m);
    for n=1:1:size(Hn,2)
        H = Hn(n);
        %solve \lambda and Ci
        fun = @(x) x*sin(x)/W-hb*cos(x)/k;
        for i=1:1:N
            x0 = (i-1)*pi;
            x1 = 0.5*pi+(i-1)*pi;
            lambda(i) = fzero(fun,[x0,x1])/W;
            C(i) = qdd*sin(lambda(i)*c)/(lambda(i)*k)/((lambda(i)*sinh(lambda(i)*H)+hbeff*cosh(lambda(i)*H)/k)*(W/2+sin(2*lambda(i)*W)/(4*lambda(i))));
        end
        for i=1:1:N
            maxt(m,n) = maxt(m,n)+C(i)*cosh(lambda(i)*H)+hbeff*C(i)*sinh(lambda(i)*H)/(k*lambda(i));
        end
    end
end
maxt = maxt+20;

figure(3)
label = {};
for i = 1:1:size(Wm,2)
    label = [label,strcat('W=',num2str(100*Wm(i)),'cm')];
    plot(Hn,maxt(i,:),'LineWidth',2)
    hold on;
end
grid on;
xlim([0 max(Hn)])
ylim([50 200])
xlabel({'$H[m]$'},'FontSize',20,'Interpreter','Latex');
ylabel({'$t_{max}[{}^{\circ}C]$'},'FontSize',20,'Interpreter','Latex');
set(gca, 'FontName','Times New Roman','FontSize', 20)
legend(label)
