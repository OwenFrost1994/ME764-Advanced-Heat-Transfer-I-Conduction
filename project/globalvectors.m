%====================== Assemble the global vectors =================
%
function [g,q,Ht] = globalvectors(nDof,nEle,nNodes,nNoEl,connArray,edgeArray,coorNo,thermal_SD,thermal_BC)
%nDof            number of degree of freedom
%nEle            number of elememts
%nNodes          number of nodes
%nNoEl           number of nodes in an elements
%connArray       connection array
%edgeArray       edge array, the node and boundary number of all edges
%coorNo          coordinates of nodes
%thermal_SD      thermal property for each domain
%thermal_BC      thermal property for each boundary
%
%Generate the element vectors(g-generation,q-heat flux,Ht-Ht_inf)
%Assemble the global vectors
%
g = zeros(nNodes,1);%empty generation vector
q = zeros(nNodes,1);%empty heat flux vector
Ht = zeros(nNodes,1);%empty Ht vector

%% generate the generation vector
lmncoord = zeros(nNoEl,nDof);
%
%Loop over all the elements to generate generation vector for each node
for lmn = 1:1:nEle
    %
    %Extract coords of nodes, DOF for the current element
    %
    for a = 1:1:nNoEl%loop of element nodes
        for i = 1:1:nDof%pick out node coordinates
          lmncoord(a,i)=coorNo(connArray(lmn,a),i);
        end
    end
    %generation vector for a single element
    gel=elgenvector(nDof,nNoEl,lmncoord,thermal_SD(connArray(lmn,nNoEl+1),2));
    
    % Add the current element generation vector to the global generation
    % vector
    for a = 1:1:nNoEl
        rw = connArray(lmn,a);
        g(rw,1) = g(rw,1) + gel(a,1);
    end
end

%% generate the heat flux vector
lmncoord = zeros(2,nDof);
%
%Loop over all the elements to generate every heat flux vector
for lmn = 1:1:size(edgeArray,1)
    %
    %Extract coords of nodes, DOF for the current element
    %
    for a = 1:1:2%loop of element nodes
        for i = 1:1:nDof%pick out node coordinates
            lmncoord(a,i)=coorNo(edgeArray(lmn,a),i);
        end
    end
    %heat flux vector for a single edge
    qel=elheatf(lmncoord,thermal_BC(edgeArray(lmn,3),3));
    
    % Add the current element capacity to the global heat flux vector
    for a = 1:1:2
        rw = edgeArray(lmn,a);
        q(rw,1) = q(rw,1) + qel(a,1);
    end
end

%% generate the Ht vector(boundary temperature)
lmncoord = zeros(2,nDof);
%
%Loop over all the elements to generate every edge vector
for lmn = 1:1:size(edgeArray,1)
    %
    %Extract coords of nodes, DOF for the current element
    %
    for a = 1:1:2%loop of element nodes
        for i = 1:1:nDof%pick out node coordinates
            lmncoord(a,i)=coorNo(edgeArray(lmn,a),i);
        end
    end
    %temperature vector for a single edge
    Htel=eltinf(lmncoord,thermal_BC(edgeArray(lmn,3),1),thermal_BC(edgeArray(lmn,3),2));
    
    % Add the current element capacity to the global heat flux vector
    for a = 1:1:2
        rw = edgeArray(lmn,a);
        Ht(rw,1) = Ht(rw,1) + Htel(a,1);
    end
end
end