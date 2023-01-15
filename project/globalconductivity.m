%====================== Assemble the global conductivity matrix =================
%
function [K] = globalconductivity(nDof,nEle,nNodes,nNoEl,connArray,coorNo,thermal_SD)
%nDof            number of degree of freedom
%nEle            number of elememts
%nNodes          number of nodes
%nNoEl           number of nodes in an elements
%connArray       connection array
%thermal_SD      thermal property for each domain
%
%Generate the element conductivity matrix
%Assemble the global conductivity matrix
%
K = zeros(nNodes,nNodes);%empty global conductivity
lmncoord = zeros(nNoEl,nDof);
%
%Loop over all the elements to generate every conductivity
for lmn = 1:nEle
%
%Extract coords of nodes, DOF for the current element
%
    for a = 1:nNoEl%loop of element nodes
        for i = 1:nDof%pick out node coordinates
          lmncoord(a,i)=coorNo(connArray(lmn,a),i);
        end
    end
    
    if connArray(lmn,nNoEl+1) == 1%element in matrix
        kel=elcond(nDof,nNoEl,lmncoord,thermal_SD(1,1));
    end
    if connArray(lmn,nNoEl+1) == 2%element in layer
        kel=elcond(nDof,nNoEl,lmncoord,thermal_SD(1,2));
    end
     
    % Add the current element conductivity to the global conductivity
    for a = 1:nNoEl
        for b = 1:nNoEl
            rw = connArray(lmn,a);
            cl = connArray(lmn,b);
            K(rw,cl) = K(rw,cl) + kel(a,b);
        end
    end
end
end