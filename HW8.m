clear
clc
close all

%% d)
theta = 0.01:0.01:10;
x = 0.01:0.01:5;

beta = [0.1,1,10];

alpha = 0;
tol = 1E-9;

figure(1)
for i=1:1:size(beta,2)
    for j=1:1:size(x,2)
        t(j,:) = invlap('that', theta', alpha, tol, x(j), beta(i));
    end
    surf(x,theta,t','EdgeColor','none')
    [row, col] = find(max(max(t))==t);
    text(x(row),theta(col),max(max(t)),strcat({'${\beta}=$'},num2str(beta(i))),'Color','red','FontSize',20,'interpreter','latex');
    hold on;
end
xlim([0 max(x)])
ylim([0 max(theta)])
set(gcf,'position',[200,300,1000,800])
xlabel({'$\bar{x}$'},'FontSize',20,'Interpreter','Latex');
ylabel({'$\bar{\theta}$'},'FontSize',20,'Interpreter','Latex');
zlabel({'$\bar{t}$'},'FontSize',20,'Interpreter','Latex');
set(gca, 'FontName','Times New Roman','FontSize', 20)