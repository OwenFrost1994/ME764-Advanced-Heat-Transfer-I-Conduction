%================= EDGE HEAT FLUX VECTOR ================================
%
function qel=elheatf(coord,qs)
%
%  Assemble the Edge heat flux vector
%  Input variables
%coord(a,i)          ith coordinate of ath node
%qs                  edge heat flux
%
%  Local variables
%npoints            Numbers integration points
%xi(i,intpt)        ith local coord of integration point
%w(intpt)           weight for integration point no. intpt
%N(a)               Shape function associated with ath node on element
%qel(row,1)       element stiffness matrix, row¡ª¡ªrow index
%
%
npoints = numberofintegrationpoints(1,2);
xi = zeros(1,1);
qel = zeros(2,1);
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
     
     %Compute the element vector
%      for a=1:1:2
%          row=a;
%          qel(row,1)=qel(row,1) + qs*N(a)*w(intpt)*ds;
%      end
     qel=[qs*ds/2;qs*ds/2];
 end
end