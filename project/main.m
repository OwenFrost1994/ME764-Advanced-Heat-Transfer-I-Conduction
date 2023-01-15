clear
clc
close all
%% workpiece geometry
W=0.03;%width of workpiece
L=W;%height of workpiece
T=0.03;%thickness of workpiece
volf=0.06;%volume fracture of layer
W1=W*volf;%width of layer in workpiece

%% material properties
km = 10;
kl = 2;

cm = 1000;
cl = 1000;

rhom = 1000;
rhol = 1000;

gddd = 1E6;

ht = 100;
Tt = 30;
hb = 100;
Tb = 40;
T0 = 0;

%% shape functions and boundary conditions
% order of element (1-Linear or 2-Quadratic)
order=2;
%%nodal degree of freedom;%Number degrees of freedom per node (2 for 2D, 3 for 3D)
nDof=2;

%sd and bc
%domian 1-matrix, 2-stiff layer
thermal_SD = [km 0 rhom cm; kl gddd rhol cl;];
thermal_BC = [0 0 0; hb Tb 0; 0 0 0; ht Tt 0;];

%% Meshing parameters
nElx_M=10; %%number of elements in x in matrix
nElx_L=5; %%number of elements in x in layer
nEly=nElx_L+nElx_M;   %%number of elements in y(the division of elements in y are same in two phase)
nElz=2;   %%number of elements in z(I abandoned to realize the 3D model)

nNoEl = numNoEl(nDof,order);
%(here ndof=ncoord, but the program allows them to be different to allow extension to plate & beam elements with C^1 continuity)

%No. elements && connectivity
[nEle,nNodes,connArray,coorNoUD,edgeArray,nEleM,nEleL]=meshgeneration(nElx_M,nElx_L,nEly,nElz,W,W1,L,T,order,nDof);
%nNodes       Number of nodes
%nEle       Number of elements
%connArray       The node number of ith element's jth node in global node marking
%coorNoUD       The coordinates of ith node in Undeformed status, jth column is the jth coordinate
%nEleM       Number of elements in matrix
%nEleL       Number of elements in layer

% Plot the initial mesh for check
meshplotFunc(nDof,nEle,nNoEl,connArray,coorNoUD,1)

%%thermal matrices
C = globalcapacity(nDof,nEle,nNodes,nNoEl,connArray,coorNoUD,thermal_SD);
K = globalconductivity(nDof,nEle,nNodes,nNoEl,connArray,coorNoUD,thermal_SD);
H = globalconvection(nDof,size(edgeArray,1),nNodes,edgeArray,coorNoUD,thermal_BC);
[g,q,Ht] = globalvectors(nDof,nEle,nNodes,nNoEl,connArray,edgeArray,coorNoUD,thermal_SD,thermal_BC);

S = K+H;
r = q+g+Ht;

tss = pinv(S)*r;
X=coorNoUD(:,1);
Y=coorNoUD(:,2);
[x,y]=meshgrid(0:W/20:W,0:L/20:L);
Tss=griddata(X,Y,tss,x,y,'v4');

figure(2)
contourf(x,y,Tss,'LineWidth',1)
colorbar
set(gcf,'position',[200,300,1000,800])
xlabel({'$x[m]$'},'FontSize',20,'Interpreter','Latex');
ylabel({'$t[m]$'},'FontSize',20,'Interpreter','Latex');
set(gca, 'FontName','Times New Roman','FontSize', 20)
