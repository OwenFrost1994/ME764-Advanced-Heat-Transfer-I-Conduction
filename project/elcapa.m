%================= ELEMENT CAPACITY MATRIX ================================
%
function [cel]=elcapa(nDof,nNoEl,coord,rho,c)
%
%  Assemble the element capacity
%  Input variables
%nDof                Numbers of degrees of freedom per node (2 in 2D and 3 in 3D)
%nNoEl               Numbers nodes on the element
%coord(a,i)          ith coordinate of ath node
%c                   Material heat cacpacity
%
%  Local variables
%npoints            Numbers integration points
%xi(i,intpt)        ith local coord of integration point
%w(intpt)           weight for integration point no. intpt
%N(a)               Shape function associated with ath node on element
%cel(row,col)       element capacity matrix, row！！row index, col！！column
%index
%
%
npoints = numberofintegrationpointscapactiy(nDof,nNoEl);
dxdxi = zeros(nDof,nDof);
xi = zeros(nDof,1);
cel = zeros(nNoEl,nNoEl);
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
     
     %Compute the element capacity
%      for a = 1:1:nNoEl
%         for b = 1:1:nNoEl
%            for i = 1:1:nDof
%               row = nDof*(a-1)+i;
%               col = nDof*(b-1)+i;
%               cel(col,row) = cel(col,row) + c*rho*N(b)*N(a)*w(intpt)*dt;
%            end
%         end
%       end
%    end
     for a=1:1:nNoEl
         for b=1:1:nNoEl
             row=a;
             col=b;
             cel(row,col)=cel(row,col) + rho*c*N(b)*N(a)*w(intpt)*dt;
         end
     end
 end
end