clear
clc
close all

%% d)
N = 10;
y = 0:N*pi/(100*N):N*pi;
for i=1:1:size(y,2)
    res(i) = besselj(0,y(i));
end

%%solve \lambda
fun = @(x) besselj(0,x);
label = {};
for i = 1:1:N
    label = [label,strcat(num2str(i),'\pi')];
    x0 = (i-0.5)*pi;
    x1 = i*pi;
    lambda(i) = fzero(fun,[x0,x1]);
end
label = ['0',label];

figure(1)
subplot(1,2,1)
plot(y,res,'b-','LineWidth',2);
hold on;
for i=1:1:N
    plot([(i-0.5)*pi,(i-0.5)*pi],[min(res),max(res)],'--g','LineWidth',2)
    hold on;
    plot([i*pi,i*pi],[min(res),max(res)],'--y','LineWidth',2)
    hold on;
end
plot(lambda,zeros(size(lambda)),'.r','Markersize',30,'LineWidth',2);
hold on;
grid on;
xlim([0 N*pi])
set(gca,'xtick',0:pi:10*pi);
set(gca,'xticklabels',label);
xlabel({'$\lambda_{i}$'},'FontSize',20,'Interpreter','Latex');
ylabel({'$res$'},'FontSize',20,'Interpreter','Latex');
set(gca, 'FontName','Times New Roman','FontSize', 20)
subplot(1,2,2)
plot(y,abs(res),'b-','LineWidth',2);
hold on;
for i=1:1:N
    plot([(i-0.5)*pi,(i-0.5)*pi],[min(res),max(res)],'--g','LineWidth',2)
    hold on;
    plot([i*pi,i*pi],[min(res),max(res)],'--y','LineWidth',2)
    hold on;
end
plot(lambda,zeros(size(lambda)),'.r','Markersize',30,'LineWidth',2);
hold on;
grid on;
xlim([0 N*pi])
set(gca,'xtick',0:pi:10*pi);
set(gca,'xticklabels',label);
xlabel({'$\lambda_{i}$'},'FontSize',20,'Interpreter','Latex');
ylabel({'$abs(res)$'},'FontSize',20,'Interpreter','Latex');
set(gca, 'FontName','Times New Roman','FontSize', 20)

N = 50;
%%solve \lambda
fun = @(x) besselj(0,x);
for i = 1:1:N
    x0 = (i-0.5)*pi;
    x1 = i*pi;
    lambda(i) = fzero(fun,[x0,x1]);
    fun1 = @(x) (x-x.^3).*besselj(0,lambda(i)*x)/4;
    fun2 = @(x) x.*besselj(0,lambda(i)*x).^2;
    C(i) = integral(fun1,0,1)/integral(fun2,0,1);
end

r = 0:1/100:1;
theta = [0,0.01,0.1,0.2,0.5,1];
v = zeros(size(theta,2),size(r,2));
for k=1:1:size(theta,2)
    for j=1:1:size(r,2)
        for i=1:1:N
            v(k,j) = v(k,j) + C(i)*besselj(0,lambda(i)*r(j))*exp(-lambda(i)^2*theta(k));
        end
    end
end

figure(2)
label = {};
for i = 1:1:size(theta,2)
    label = [label,strcat('$\bar{\theta}=$',num2str(theta(i)))];
    plot(r,v(i,:),'LineWidth',2)
    hold on;
end
grid on;
xlim([0 max(r)])
xlabel({'$\bar{r}$'},'FontSize',20,'Interpreter','Latex');
ylabel({'$v(\bar{r},\bar{\theta})$'},'FontSize',20,'Interpreter','Latex');
set(gca, 'FontName','Times New Roman','FontSize', 20)
legend(label,'Interpreter','Latex')

%% e)
U = zeros(size(theta,2),size(r,2));
for k=1:1:size(theta,2)
    for j=1:1:size(r,2)
        U(k,j) = U(k,j) + (1-r(j)^2)/4;
        for i=1:1:N
            U(k,j) = U(k,j) - C(i)*besselj(0,lambda(i)*r(j))*exp(-lambda(i)^2*theta(k));
        end
    end
end

figure(3)
label = {};
for i = 1:1:size(theta,2)
    label = [label,strcat('$\bar{\theta}=$',num2str(theta(i)))];
    plot(r,U(i,:),'LineWidth',2)
    hold on;
end
grid on;
xlim([0 max(r)])
xlabel({'$\bar{r}$'},'FontSize',20,'Interpreter','Latex');
ylabel({'$U(\bar{r},\bar{\theta})$'},'FontSize',20,'Interpreter','Latex');
set(gca, 'FontName','Times New Roman','FontSize', 20)
legend(label,'Interpreter','Latex')

%% f)
theta = [0,0.1,0.2,0.5,1,2,5];
t = zeros(size(theta,2),size(r,2));
beta = 1;
for k=1:1:size(theta,2)
    for j=1:1:size(r,2)
        t(k,j) = t(k,j) + exp(-beta*theta(k))*(1-r(j)^2)/4;
        for i=1:1:N
            t(k,j) = t(k,j) - beta*C(i)*besselj(0,lambda(i)*r(j))*(exp(-beta*theta(k))-exp(-lambda(i)^2*theta(k)))/(lambda(i)^2-beta);
        end
    end
end

figure(4)
label = {};
for i = 1:1:size(theta,2)
    label = [label,strcat('$\bar{\theta}=$',num2str(theta(i)))];
    plot(r,t(i,:),'LineWidth',2)
    hold on;
end
grid on;
xlim([0 max(r)])
xlabel({'$\bar{r}$'},'FontSize',20,'Interpreter','Latex');
ylabel({'$\bar{t}(\bar{r},\bar{\theta})$'},'FontSize',20,'Interpreter','Latex');
set(gca, 'FontName','Times New Roman','FontSize', 20)
legend(label,'Interpreter','Latex')

beta = [0.1,0.2,0.5,1,1.5,2,3,4,5];
theta = [0,0.1,0.2,0.5,1,2];
maxt = zeros(size(theta,2),size(beta,2));
for k=1:1:size(beta,2)
    for j=1:1:size(theta,2)
        maxt(j,k) = maxt(j,k) + exp(-beta(k)*theta(j))*(1-0^2)/4;
        for i=1:1:N
            maxt(j,k) = maxt(j,k) - beta(k)*C(i)*besselj(0,lambda(i)*0)*(exp(-beta(k)*theta(j))-exp(-lambda(i)^2*theta(j)))/(lambda(i)^2-beta(k));
        end
    end
end

figure(5)
for i = 1:1:size(theta,2)
    label = [label,strcat('$\bar{\theta}=$',num2str(theta(i)))];
    plot(beta,maxt(i,:),'LineWidth',2)
    hold on;
end
grid on;
xlim([0 max(beta)])
xlabel({'$\beta$'},'FontSize',20,'Interpreter','Latex');
ylabel({'$\bar{t}_{max}$'},'FontSize',20,'Interpreter','Latex');
set(gca, 'FontName','Times New Roman','FontSize', 20)
legend(label,'Interpreter','Latex')
title('Maximum temperature under different \beta')