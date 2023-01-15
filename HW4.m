clear
clc
close all

th = 0.01;
Rw = 0.03;
hb = 50;
Tedge = 25+273;
Tinf = 20+273;
k = 1.2;
qrad = 1000;
alpha = 100;

%% a)
% the intervals of eigenvalue \lambda and the values of \lambda
N = 10;
x = 0:N*pi/(100*N):N*pi;
for i = 1:1:size(x,2)
    res(i) = x(i)*sin(x(i))/th-hb*cos(x(i))/k;
end
label = {};
for i = 1:1:N
    label = [label,strcat(num2str(i),'\pi')];
end
label = ['0',label];

%%solve \lambda
fun = @(x) x*sin(x)-hb*th*cos(x)/k;
for i=1:1:N
    x0 = (i-1)*pi;
    x1 = 0.5*pi+(i-1)*pi;
    lambda(i) = fzero(fun,[x0,x1])/th;
end

figure(1)
subplot(1,2,1)
plot(x,res,'b-','LineWidth',2);
hold on;
for i = 1:1:N-1
    plot([0.5*pi+(i-1)*pi,0.5*pi+(i-1)*pi],[min(res),max(res)],'--k','LineWidth',2)
    hold on;
end
plot(lambda*th,zeros(size(lambda)),'.r','Markersize',30,'LineWidth',2);
hold on;
grid on;
xlim([0 N*pi])
set(gca,'xtick',0:pi:10*pi);
set(gca,'xticklabels',label);
xlabel({'$\lambda_{i}th$'},'FontSize',20,'Interpreter','Latex');
ylabel({'$res$'},'FontSize',20,'Interpreter','Latex');
set(gca, 'FontName','Times New Roman','FontSize', 20)
subplot(1,2,2)
plot(x,abs(res),'b-','LineWidth',2);
hold on;
for i = 1:1:N-1
    plot([0.5*pi+(i-1)*pi,0.5*pi+(i-1)*pi],[min(abs(res)),max(abs(res))],'--k','LineWidth',2)
    hold on;
end
plot(lambda*th,zeros(size(lambda)),'.r','Markersize',30,'LineWidth',2);
hold on;
grid on;
xlim([0 N*pi])
set(gca,'xtick',0:pi:10*pi);
set(gca,'xticklabels',label);
xlabel({'$\lambda_{i}th$'},'FontSize',20,'Interpreter','Latex');
ylabel({'$abs(res)$'},'FontSize',20,'Interpreter','Latex');
set(gca, 'FontName','Times New Roman','FontSize', 20)

%% b)
% C4 for a(x)
C4 = Tinf + qrad*(exp(-alpha*th)+alpha*th)/(alpha*k)-qrad*(exp(-alpha*th)-1)/hb;
N = 50;
%%solve \lambda
fun = @(x) x*sin(x)-hb*th*cos(x)/k;
for i=1:1:N
    x0 = (i-1)*pi;
    x1 = 0.5*pi+(i-1)*pi;
    lambda(i) = fzero(fun,[x0,x1])/th;
end

%%calculate Ai
for i=1:1:N
    x0 = 0;
    x1 = th;
    fun1 = @(x) (Tedge-(-qrad*(exp(-alpha*x)+alpha*x)/(alpha*k)+C4)).*cos(lambda(i)*x);
    A(i) = integral(fun1,x0,x1)/(besseli(0,lambda(i)*Rw)*(th/2+sin(2*lambda(i)*th)/(4*lambda(i))));
end

x = 0:th/5:th;
r = 0:Rw/100:Rw;
t = zeros(size(x,2),size(r,2));
for i=1:1:size(x,2)
    a = -qrad*(exp(-alpha*x(i))+alpha*x(i))/(alpha*k)+C4;
    for j=1:1:size(r,2)
        for m=1:1:N
            t(i,j) = t(i,j)+A(m)*cos(lambda(m)*x(i))*besseli(0,lambda(m)*r(j));
        end
        t(i,j) = t(i,j)+a-273;
    end
end
%
figure(2)
for i = 1:1:6
    plot(r,t(i,:),'LineWidth',2)
    hold on;
end
grid on;
xlim([0 Rw])
xlabel({'$r$'},'FontSize',20,'Interpreter','Latex');
ylabel({'$t[{}^{\circ}C]$'},'FontSize',20,'Interpreter','Latex');
set(gca, 'FontName','Times New Roman','FontSize', 20)
legend('x=0cm','x=0.2cm','x=0.4cm','x=0.6cm','x=0.8cm','x=1cm')


x = 0:th/100:th;
r = 0:Rw/1000:Rw;
t = zeros(size(x,2),size(r,2));
for i=1:1:size(x,2)
    a = -qrad*(exp(-alpha*x(i))+alpha*x(i))/(alpha*k)+C4;
    for j=1:1:size(r,2)
        for m=1:1:N
            t(i,j) = t(i,j)+A(m)*cos(lambda(m)*x(i))*besseli(0,lambda(m)*r(j));
        end
        t(i,j) = t(i,j)+a-273;
    end
end
figure(3)
contour(r,x,t)
xlabel({'$r$'},'FontSize',20,'Interpreter','Latex');
ylabel({'$x$'},'FontSize',20,'Interpreter','Latex');
set(gca, 'FontName','Times New Roman','FontSize', 20)
title('Temperature distribution in the window')