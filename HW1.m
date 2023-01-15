%% c)
a = 1;
b = 2;

x1 = 0:0.01:a;
x2 = a:0.01:a+b;

t1 = (tanh(a+b)*cosh(a)-sinh(a))*sinh(x1);
t2 = cosh(a)*(tanh(a+b)*sinh(x2)-cosh(x2))+1;

figure(1)
plot([x1,x2],[t1,t2],'b-','Markersize',6,'LineWidth',2)
hold on;
grid on;
xlabel({'$\bar{x}$'},'FontSize',20,'Interpreter','Latex');
ylabel({'$\bar{t}$'},'FontSize',20,'Interpreter','Latex');
set(gca, 'FontName','Times New Roman','FontSize', 20)
legend(strcat('a=',num2str(a),',b=',num2str(b)));

%% d)
b = [0:0.2:1, 2:1:5];
a = 0:0.01:5;

figure(2)
for j = 1:1:size(b,2)
    for i = 1:1:size(a,2)
        tmax(j,i) = cosh(a(i))*(tanh(a(i)+b(j))*sinh(a(i)+b(j))-cosh(a(i)+b(j)))+1;
    end
    plot(a,tmax(j,:),'LineWidth',2);
    text(max(a(floor(size(a,2)/2))), max(tmax(j,:)), strcat('b=',num2str(b(j))),'FontSize',20);
    hold on;
end
grid on;
xlabel({'$a$'},'FontSize',20,'Interpreter','Latex');
ylabel({'$\bar{t}_{max}$'},'FontSize',20,'Interpreter','Latex');
set(gca, 'FontName','Times New Roman','FontSize', 20)

%% e)
th = 0.001;
hb = 20;
k = 10;
L1 = 0.1;
L2 = 0.05;
q = 1000;
tinf = 20;

L = sqrt(k*th/(2*hb));
T = q/hb;

a = L1/L;
b = L2/L;

x1 = 0:0.01:a;
x2 = a:0.01:a+b;

t1 = (tanh(a+b)*cosh(a)-sinh(a))*sinh(x1);
t2 = cosh(a)*(tanh(a+b)*sinh(x2)-cosh(x2))+1;

x = L*[x1,x2];
t = [t1,t2]*T + tinf;

figure(3)
plot(x,t,'LineWidth',2);
hold on;
grid on;
xlabel({'$x[m]$'},'FontSize',20,'Interpreter','Latex');
ylabel({'$t[{}^{\circ}C]$'},'FontSize',20,'Interpreter','Latex');
set(gca, 'FontName','Times New Roman','FontSize', 20)