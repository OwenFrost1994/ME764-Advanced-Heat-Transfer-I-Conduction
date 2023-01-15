%====================== Assemble the global convection matrix =================
%
function [H] = globalconvection(nDof,nEdge,nNodes,edgeArray,coorNo,thermal_BC)
%nDof            number of degree of freedom
%nEdge            number of elememts
%nNodes          number of nodes
%nNoEl           number of nodes in an elements
%edgeArray       edge array, the node and boundary number of all edges
%coorNo          coordinates of nodes
%thermal_BC      thermal property for each boundary
%
%Generate the element conduction matrix
%Assemble the global conduction matrix
%
H = zeros(nNodes,nNodes);%empty global capacity
lmncoord = zeros(2,nDof);
%
%Loop over all the elements to generate every convection matrix
for lmn = 1:1:nEdge
    %
    %Extract coords of nodes, DOF for the current element
    %
    for a = 1:1:2%loop of element nodes
        for i = 1:1:nDof%pick out node coordinates
            lmncoord(a,i)=coorNo(edgeArray(lmn,a),i);
        end
    end
    %conduction matrix for a single edge
    hel=elconv(lmncoord,thermal_BC(edgeArray(lmn,3),1));
    
    % Add the current element capacity to the global capcacity matrix
    for a = 1:1:2
        for b = 1:1:2
            rw = edgeArray(lmn,a);
            cl = edgeArray(lmn,b);
            H(rw,cl) = H(rw,cl) + hel(a,b);
        end
    end
end
end