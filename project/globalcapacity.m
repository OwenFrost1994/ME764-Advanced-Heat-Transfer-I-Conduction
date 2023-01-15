%====================== Assemble the global capacity matrix =================
%
function [C] = globalcapacity(nDof,nEle,nNodes,nNoEl,connArray,coorNo,thermal_SD)
%nDof            number of degree of freedom
%nEle            number of elememts
%nNodes          number of nodes
%nNoEl           number of nodes in an elements
%connArray       connection array
%thermal_SD      thermal property for each domain
%
%Generate the element capacity matrix
%Assemble the global capacity matrix
%
C = zeros(nNodes,nNodes);%empty global capacity
lmncoord = zeros(nNoEl,nDof);
%
%Loop over all the elements to generate every conductivity
for lmn = 1:1:nEle
    %
    %Extract coords of nodes, DOF for the current element
    %
    for a = 1:1:nNoEl%loop of element nodes
        for i = 1:1:nDof%pick out node coordinates
          lmncoord(a,i)=coorNo(connArray(lmn,a),i);
        end
    end
    %capacity matrix for a single element
    if connArray(lmn,nNoEl+1) == 1%element in matrix
        cel=elcapa(nDof,nNoEl,lmncoord,thermal_SD(1,3),thermal_SD(1,4));
%         cel=elcapa(nDof,nNoEl,coord,rho,c);
    end
    if connArray(lmn,nNoEl+1) == 2%element in layer
        cel=elcapa(nDof,nNoEl,lmncoord,thermal_SD(2,3),thermal_SD(2,4));
%         cel=elcapa(nDof,nNoEl,coord,rho,c);
    end
    
    % Add the current element capacity to the global capcacity matrix
    for a = 1:1:nNoEl
        for b = 1:1:nNoEl
            rw = connArray(lmn,a);
            cl = connArray(lmn,b);
            C(rw,cl) = C(rw,cl) + cel(a,b);
        end
    end
end
end