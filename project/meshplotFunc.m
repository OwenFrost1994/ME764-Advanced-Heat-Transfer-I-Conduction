function meshplotFunc(nDof,nEle,nNoEl,connArray,coorNoUD,plotn)
%nDof    node degree of freedom 2-2D,3-3D
%nEle    number of elememts
%nNodes    number of nodes
%nNoEl    number of nodes in an elements
%connArray    connection array
%coorNoUD    coordinates of ith node in undeformed status
%order    order of element
%W,H,W1    width, hight, width of layer
%plotn    number of plot
% Function to plot a mesh.  
f2D_4 = [1,2,3,4];
f2D_8 = [1,5,2,6,3,7,4,8];
f3D_8 = [[1,2,3,4];[5,8,7,6];[1,5,6,2];[2,3,7,6];[3,7,8,4];[4,8,5,1]];
f3D_20 = [[1,9,2,10,3,11,4,12];[5,16,8,15,7,14,6,13];
    [1,17,5,13,6,18,2,9];[2,18,6,14,7,19,3,10];
    [3,19,7,15,8,20,4,11];[4,20,8,16,5,17,1,12]];
figure(plotn)
if (nDof==2)  % Plot a 2D mesh
    for lmn=1:1:nEle
        for i=1:1:nNoEl
            x(i,1:2) = coorNoUD(connArray(lmn,i),1:2);
        end
        scatter(x(:,1),x(:,2),24,'r','LineWidth',2.0);
        hold on;
        if nNoEl==4
            patch('Vertices',x,'Faces',f2D_4,'FaceColor','none','EdgeColor','b','LineWidth',2.0);
            hold on;
        else if nNoEl==8
                patch('Vertices',x,'Faces',f2D_8,'FaceColor','none','EdgeColor','b','LineWidth',2.0);
                hold on;
            end
        end
    end
    grid on;
    xlabel('x(m)')
    ylabel('y(m)')
    set(gca, 'FontName','Times New Roman','FontSize', 20)
    set(get(gca,'XLabel'),'Fontsize',20)
    set(get(gca,'YLabel'),'Fontsize',20)
    title('Undeformed Mesh');
else% Plot a 3D mesh%ÉÐÎ´debug
    for lmn=1:1:nEle
        for i=1:1:nNoEl
            x(i,1:3) = coorNoUD(connArray(lmn,i),1:3);
        end
        scatter3(x(:,1),x(:,2),x(:,3),24,'r','LineWidth',2.0);
        if (nNoEl==8) 
            patch('Vertices',x,'Faces',f3D_8,'FaceColor','none','EdgeColor','b','LineWidth',2.0);
            hold on;
        else if (nNoEl==20)
                patch('Vertices',x,'Faces',f3D_20,'FaceColor','none','EdgeColor','b','LineWidth',2.0);
                hold on;
            end
        end
    end
    grid on;
    xlabel('x(m)')
    ylabel('y(m)')
    zlabel('z(m)')
    view([0.5,0.5,0.5])
    set(gca, 'FontName','Times New Roman','FontSize', 20)
    set(get(gca,'XLabel'),'Fontsize',20)
    set(get(gca,'YLabel'),'Fontsize',20)
    title('Mesh');
end
axis equal
hold off
end