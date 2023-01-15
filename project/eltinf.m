%================= EDGE TEMPERATURE VECTOR ================================
%
function Htel=eltinf(coord,h,Tinf)
%
%  Assemble the Edge heat flux vector
%  Input variables
%coord(a,i)          ith coordinate of ath node
%h                   boundary convection
%Tinf                boundary temperature
%
%  Local variables
%npoints            Numbers integration points
%xi(i,intpt)        ith local coord of integration point
%w(intpt)           weight for integration point no. intpt
%N(a)               Shape function associated with ath node on element
%Htel(row,1)       element stiffness matrix, row¡ª¡ªrow index
%
%
npoints = numberofintegrationpoints(1,2);
xi = zeros(1,1);
Htel = zeros(2,1);
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
%          Htel(row,1)=Htel(row,1) + h*Tinf*N(a)*w(intpt)*ds;
%      end
     Htel=[h*Tinf*ds/2;h*Tinf*ds/2];
 end
end