%================= ELEMENT CONVECTION MATRIX ================================
%
function [cel]=elconv(coord,h)
%
%  Assemble the element capacity
%  Input variables
%nDof                Numbers of degrees of freedom per node (2 in 2D and 3 in 3D)
%nNoEl               Numbers nodes on the element
%coord(a,i)          ith coordinate of ath node
%h                   boundary convection
%
%  Local variables
%npoints            Numbers integration points
%xi(i,intpt)        ith local coord of integration point
%w(intpt)           weight for integration point no. intpt
%N(a)               Shape function associated with ath node on element
%cel(row,col)       element convection matrix, row！！row index, col！！column
%index
%
%
npoints = numberofintegrationpoints(1,2);
xi = zeros(1,1);
cel = zeros(2,2);
%
%  Set up integration points && weights    
%
xilist=integrationpoints(1,2,npoints);
w=integrationweights(1,2,npoints);
%
%  Loop over the integration points(intpt)
%
 for intpt=1:1:npoints
     %Compute shape functions && derivatives wrt local coords
     for i=1:1:1
         xi(i) = xilist(i,intpt);
     end
     N = shapefunctions(2,1,xi);
     ds = sqrt((coord(1,1)-coord(2,1))^2+(coord(1,2)-coord(2,2))^2);
     
     %Compute the element capacity
%      for a=1:1:2
%          for b=1:1:2
%              row=a;
%              col=b;
%              cel(row,col)=cel(row,col) + h*N(b)*N(a)*w(intpt)*ds;
%          end
%      end
     cel=[h*ds/3,h*ds/6;h*ds/6,h*ds/3];
 end
end