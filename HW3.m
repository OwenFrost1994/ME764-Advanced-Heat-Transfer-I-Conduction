close all
clc
clear all
%% f)
beta = 1;
N = 100;
xb = 0:0.01:1;

thetab = [0.0001,0.001,0.01,0.1,0.2,0.5,1,10,100];
a = -(sinh(beta*xb)-tanh(beta)*cosh(beta*xb))/beta;

figure(1)
for j=1:1:size(thetab,2)
    u = 0;
    for i=1:1:N
        lambda = (2*i-1)*pi/2;
        u = u - 2*cos(lambda*xb)*exp(-(beta^2+lambda^2)*thetab(j))/(beta^2+lambda^2);
    end
    t = u+a;
    plot(xb,t,'LineWidth',2);
    text(xb(1), max(t), strcat({'$\bar{\theta}=$'},num2str(thetab(j))),'FontSize',20,'interpreter','latex');
    hold on;
end
grid on;
xlabel({'$\bar{x}$'},'FontSize',20,'Interpreter','Latex');
ylabel({'$\bar{t}$'},'FontSize',20,'Interpreter','Latex');
set(gca, 'FontName','Times New Roman','FontSize', 20)

%% g)
beta = [0.1:0.2:1,1:1:10];
N = [1,5,10:10:100];
xb = 0;

thetabm = [0.01,0.1,1,10];
for m=1:1:size(thetabm,2)
    thetab = thetabm(m);
    figure(1+m)
    for k=1:1:size(N,2)
        u = zeros(size(beta));
        a = zeros(size(beta));
        for j=1:1:size(beta,2)
            a(j) = -(sinh(beta(j)*xb)-tanh(beta(j))*cosh(beta(j)*xb))/beta(j);
            for i=1:1:N(k)
                lambda = (2*i-1)*pi/2;
                u(j) = u(j) - 2*cos(lambda*xb)*exp(-(beta(j)^2+lambda^2)*thetab)/(beta(j)^2+lambda^2);
            end
        end
        t = u+a;
        plot(beta,t,'LineWidth',2);
        text(beta(1), max(t), strcat({'${N}=$'},num2str(N(k))),'FontSize',20,'interpreter','latex');
        hold on;
    end
    grid on;
    xlabel({'${\beta}$'},'FontSize',20,'Interpreter','Latex');
    ylabel({'$\bar{t}_{max}$'},'FontSize',20,'Interpreter','Latex');
    set(gca, 'FontName','Times New Roman','FontSize', 20)
    title(strcat({'$\bar{\theta}=$'},num2str(thetab)),'FontSize',20,'interpreter','latex')
end
%% 
beta = [0.01:0.01:0.1,0.2:0.1:1,2:1:10,20:10:100,200:100:1000];
lambda = pi/2;
for i=1:1:size(beta,2)
    tbmax(i) = tanh(beta(i))/beta(i);
    theta(i) = -log(tanh(beta(i))*(lambda^2+beta(i)^2)/(20*beta(i)))/(lambda^2+beta(i)^2);
end
figure(2+size(thetabm,2))
plot(beta,tbmax,'r-','LineWidth',2);
hold on;
plot(beta,theta,'b-','LineWidth',2);
hold on;
grid on;
set(gca,'xscale','log')
xlabel({'${\beta}$'},'FontSize',20,'Interpreter','Latex');
ylabel({'$\bar{t}_{max}$'},'FontSize',20,'Interpreter','Latex');
set(gca, 'FontName','Times New Roman','FontSize', 20)
legend({'$\bar{t}_{max}$','$\bar{\theta}_{90\%\bar{t}_{max}}$'},'FontSize',20,'interpreter','latex')