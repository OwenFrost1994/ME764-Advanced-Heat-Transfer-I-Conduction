%================= ELEMENT GENERATION VECTOR ================================
%
function [gel]=elgenvector(nDof,nNoEl,coord,g)
%
%  Assemble the element capacity
%  Input variables
%nDof                Numbers of degrees of freedom per node (2 in 2D and 3 in 3D)
%nNoEl               Numbers nodes on the element
%coord(a,i)          ith coordinate of ath node
%g                   volumn generation
%
%  Local variables
%npoints            Numbers integration points
%xi(i,intpt)        ith local coord of integration point
%w(intpt)           weight for integration point no. intpt
%N(a)               Shape function associated with ath node on element
%gel(row,1)       element stiffness matrix, row¡ª¡ªrow index
%
%
npoints = numberofintegrationpoints(nDof,nNoEl);
xi = zeros(nDof,1);
gel = zeros(nNoEl,1);
%
%  Set up integration points && weights    
%
xilist=integrationpoints(nDof,nNoEl,npoints);
w=integrationweights(nDof,nNoEl,npoints);
%
%  Loop over the integration points(intpt)
%
 for intpt=1:1:npoints
     %Compute shape functions && derivatives wrt local coords
     for i=1:1:nDof
         xi(i) = xilist(i,intpt);
     end
     N = shapefunctions(nNoEl,nDof,xi);
     dNdxi = shapefunctionderivs(nNoEl,nDof,xi);
     %Compute the jacobian matrix and its determinant
     for i = 1:1:nDof
         for j = 1:1:nDof
             dxdxi(i,j) = 0.;
             for a = 1:1:nNoEl
                 dxdxi(i,j) = dxdxi(i,j) + coord(a,i)*dNdxi(a,j);
             end
         end
     end
     dt = det(dxdxi);
     
     %Compute the element vector
     for a=1:1:nNoEl
         row=a;
         gel(row,1)=gel(row,1) + g*N(a)*w(intpt)*dt;
     end
 end
end