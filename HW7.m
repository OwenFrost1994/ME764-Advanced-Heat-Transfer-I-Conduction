clear
clc
close all

W = 0.05;
H = 0.0025;
c = 0.005;
qdd = 6.5E3;
k = 1.5;
alpha = 1E-3;

%% c)
h = 250;

N = 50;
%%\lambda Ci
for i = 1:1:N
    lambda(i) = (2*i-1)*pi/(2*W);
    C(i) = qdd*sin(lambda(i)*c)/(lambda(i)*k*W)/(lambda(i)*sinh(lambda(i)*H)+h*cosh(lambda(i)*H)/k);
end

x = 0:W/100:W;
y = 0:H/20:H;
v = zeros(size(y,2),size(x,2));
for m=1:1:size(y,2)
    for j=1:1:size(x,2)
        for i=1:1:N
            v(m,j) = v(m,j) + C(i)*cos(lambda(i)*x(j))*(cosh(lambda(i)*y(m))+h*sinh(lambda(i)*y(m))/(k*lambda(i)));
        end
    end
end

figure(1)
subplot(1,2,1)
contourf(x,y,v,'LineWidth',0)
xlim([0 W])
ylim([0 H])
xlabel({'$x[m]$'},'FontSize',20,'Interpreter','Latex');
ylabel({'$y[m]$'},'FontSize',20,'Interpreter','Latex');
ylabel({'$v(x,y)$'},'FontSize',20,'Interpreter','Latex');
set(gca, 'FontName','Times New Roman','FontSize', 20)
subplot(1,2,2)
surf(x,y,v,'AlphaData',0)
view([1,-1,1])
xlim([0 W])
ylim([0 H])
xlabel({'$x[m]$'},'FontSize',20,'Interpreter','Latex');
ylabel({'$y[m]$'},'FontSize',20,'Interpreter','Latex');
zlabel({'$v(x,y)$'},'FontSize',20,'Interpreter','Latex');
set(gca, 'FontName','Times New Roman','FontSize', 20)

%% h)
h = 250;
w = 10;

%%\lambda
for i = 1:1:N
    lambda(i) = (2*i-1)*pi/(2*W);
    lambda1(i) = sqrt(lambda(i)^2-1i*w/alpha);
    A(i) = qdd*1i*sin(lambda(i)*c)/(lambda(i)*k*W)/(lambda1(i)*sinh(lambda1(i)*H)+h*cosh(lambda1(i)*H)/k);
end

x = 0:W/100:W;
y = 0:H/20:H;
B = zeros(size(y,2),size(x,2));
for m=1:1:size(y,2)
    for j=1:1:size(x,2)
        for i=1:1:N
            B(m,j) = B(m,j) + A(i)*cos(lambda(i)*x(j))*(cosh(lambda1(i)*y(m))+h*sinh(lambda1(i)*y(m))/(k*lambda1(i)));
        end
    end
end
figure(2)
subplot(1,2,1)
contourf(x,y,abs(B),'LineWidth',0)
xlim([0 W])
ylim([0 H])
xlabel({'$x[m]$'},'FontSize',20,'Interpreter','Latex');
ylabel({'$y[m]$'},'FontSize',20,'Interpreter','Latex');
ylabel({'$|u_{cc}|$'},'FontSize',20,'Interpreter','Latex');
set(gca, 'FontName','Times New Roman','FontSize', 20)
subplot(1,2,2)
surf(x,y,abs(B),'AlphaData',0)
view([1,-1,1])
xlim([0 W])
ylim([0 H])
xlabel({'$x[m]$'},'FontSize',20,'Interpreter','Latex');
ylabel({'$y[m]$'},'FontSize',20,'Interpreter','Latex');
zlabel({'$|u_{cc}|$'},'FontSize',20,'Interpreter','Latex');
set(gca, 'FontName','Times New Roman','FontSize', 20)

%% i)
h = 250;
w = [10,100,1000];

x = 0:W/100:W;
y = 0:H/20:H;

figure(3)
for n=1:1:size(w,2)
    %%\lambda
    for i = 1:1:N
        lambda(i) = (2*i-1)*pi/(2*W);
        lambda1(i) = sqrt(lambda(i)^2-1i*w(n)/alpha);
        A(i) = qdd*1i*sin(lambda(i)*c)/(lambda(i)*k*W)/(lambda1(i)*sinh(lambda1(i)*H)+h*cosh(lambda1(i)*H)/k);
    end
    B = zeros(size(y,2),size(x,2));
    for m=1:1:size(y,2)
        for j=1:1:size(x,2)
            for i=1:1:N
                B(m,j) = B(m,j) + A(i)*cos(lambda(i)*x(j))*(cosh(lambda1(i)*y(m))+h*sinh(lambda1(i)*y(m))/(k*lambda1(i)));
            end
        end
    end
    
    surf(x,y,abs(B),'AlphaData',0)
    text(0,0,abs(B(1,1)),strcat({'${\omega}=$'},num2str(w(n))),'Color','red','FontSize',20,'interpreter','latex');
    hold on;
end
view([1,-1,1])
xlim([0 W])
ylim([0 H])
xlabel({'$x[m]$'},'FontSize',20,'Interpreter','Latex');
ylabel({'$y[m]$'},'FontSize',20,'Interpreter','Latex');
zlabel({'$|u_{cc}|$'},'FontSize',20,'Interpreter','Latex');
set(gca, 'FontName','Times New Roman','FontSize', 20)

%% i)
h = [1,5,10,15,20,50,100,200,400,800,1000];
w = [0.01:0.01:0.1,0.2:0.1:1,2:1:10,20:10:100,200:100:1000];
beta = zeros(size(h,2),size(w,2));
for m=1:1:size(h,2)
    for n=1:1:size(w,2)
        %%\lambda
    for i = 1:1:N
        lambda(i) = (2*i-1)*pi/(2*W);
        lambda1(i) = sqrt(lambda(i)^2-1i*w(n)/alpha);
        A(i) = qdd*1i*sin(lambda(i)*c)/(lambda(i)*k*W)/(lambda1(i)*sinh(lambda1(i)*H)+h(m)*cosh(lambda1(i)*H)/k);
    end
    t0 = 0;
    for i=1:1:N
        t0 = t0 + A(i)*cos(lambda(i)*0)*(cosh(lambda1(i)*0)+h(m)*sinh(lambda1(i)*0)/(k*lambda1(i)));
    end
    beta(m,n) = atan(imag(t0)/real(t0))*180/pi;%transfer into degree
    end
end

label = {};
figure(4)
for i = 1:1:size(h,2)
    label = [label,strcat('$h=$',num2str(h(i)))];
    plot(w,beta(i,:),'LineWidth',2)
    hold on;
end
grid on;
xlim([0 max(w)])
set(gcf,'position',[200,300,800,600])
set(gca,'xscale','log')
xlabel({'$\omega[rad/s]$'},'FontSize',20,'Interpreter','Latex');
ylabel({'$\beta[{}^{\circ}]$'},'FontSize',20,'Interpreter','Latex');
set(gca, 'FontName','Times New Roman','FontSize', 20)
legend(label,'Interpreter','Latex', 'Location','northwest')
title('The angle of $u_{cc}$ as function of $\omega$ and $h$','FontSize',20,'Interpreter','Latex')