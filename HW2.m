%% c)
m = [0.1:0.2:1,2:1:5,50];
r = 0:0.01:1;

figure(1)
for i =1:1:size(m,2)
    tb(i,:)=-(4+m(i)^2)*besseli(0,m(i)*r)/(m(i)^2*besseli(0,m(i)))+4/m(i)^2+r.^2;
    plot(r,tb(i,:),'LineWidth',2);
    text(m(floor(size(m,2)/2)), max(tb(i,:)), strcat('m=',num2str(m(i))),'FontSize',20);
    hold on;
end
grid on;
xlabel({'$\bar{r}$'},'FontSize',20,'Interpreter','Latex');
ylabel({'$\bar{t}$'},'FontSize',20,'Interpreter','Latex');
set(gca, 'FontName','Times New Roman','FontSize', 20)

%% d)
m = [0.1:0.1:1,2:1:10,20:10:100];
r = 0:0.01:1;

for i =1:1:size(m,2)
    tb=-(4+m(i)^2)*besseli(0,m(i)*r)/(m(i)^2*besseli(0,m(i)))+4/m(i)^2+r.^2;
    [tmax(i),rI]=max(tb);
    rmax(i)=r(rI);
end

figure(2)
plot(m,tmax,'-b','LineWidth',2);
hold on;
plot(m,rmax,'-r','LineWidth',2);
hold on;
grid on;
set(gca,'xscale','log')
xlabel({'$\bar{r}$'},'FontSize',20,'Interpreter','Latex');
legend({'$\bar{t}_{max}$','$\bar{r}_{\bar{r}_{max}}$'},'Interpreter','Latex')
set(gca, 'FontName','Times New Roman','FontSize', 20)

%% e)
th = 0.002;
r0 = 0.075;
hb = 100;
a = 1E6;
k = 1.5;
tinf = 20;

m = sqrt(hb/(k*th))*r0;
r = 0:0.01:1;

tb=-(4+m^2)*besseli(0,m*r)/(m^2*besseli(0,m))+4/m^2+r.^2;

t = a*tb*r0^2/hb+tinf;
r = r*r0;

figure(3)
plot(r,t,'-b','LineWidth',2);
hold on;
grid on;
xlabel({'$r[m]$'},'FontSize',20,'Interpreter','Latex');
ylabel({'$t[{}^{\circ}C]$'},'FontSize',20,'Interpreter','Latex');
set(gca, 'FontName','Times New Roman','FontSize', 20)