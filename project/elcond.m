%
%================= ELEMENT CONDUCTIVITY MATRIX ================================
%
function kel = elcond(nDof,nNoEl,lmncoord,k)
%
%  Assemble the element conductivity
%  Input variables
%nDof                Numbers of degrees of freedom per node (2 in 2D and 3 in 3D)
%nNoEl               Numbers nodes on the element
%lmncoord(a,i)       ith coordinate of ath node
%k                   Material heat cacpacity
%
%   Local variables
%      npoints            No. integration points
%      xi(i,inpt)         ith local coord of integration point no. intpt
%      w(intpt)           weight for integration point no. intpt
%      N(a)               Shape function associated with ath node on element
%      dNdxi(a,i)         Derivative of ath shape function wrt ith local coord
%      dNdx(a,i)          Derivative of ath shape function wrt ith global coord
%      dxdxi(i,j)         Derivative of ith global coord wrt jth local coord
%      dxidx(i,j)         Derivative of ith local coord wrt jth global coord
%      det                Determinant of jacobian
%      kel(row,col)       Rows && cols of element conductivity
%
%
npoints = numberofintegrationpoints(nDof,nNoEl);
xi = zeros(nDof,1);
dNdx = zeros(nNoEl,nDof);
dxdxi = zeros(nDof,nDof);
kel = zeros(nNoEl,nNoEl);
%
%  Set up integration points && weights    
%
xilist = integrationpoints(nDof,nNoEl,npoints);
w = integrationweights(nDof,nNoEl,npoints);
%
%  Loop over the integration points
%
for intpt = 1:1:npoints
    %Compute shape functions && derivatives wrt local coords
    for i = 1:1:nDof
        xi(i) = xilist(i,intpt);
    end
    N = shapefunctions(nNoEl,nDof,xi);
    dNdxi = shapefunctionderivs(nNoEl,nDof,xi);
    %Compute the jacobian matrix && its determinant
    for i = 1:1:nDof
        for j = 1:1:nDof
            dxdxi(i,j) = 0.;
            for a = 1:1:nNoEl
                    dxdxi(i,j) = dxdxi(i,j) + lmncoord(a,i)*dNdxi(a,j);
            end
        end
    end
    
    dxidx = inv(dxdxi);
    dt = det(dxdxi);
    %Convert shape function derivatives:derivatives wrt global coords
    for a = 1:1:nNoEl
        for i = 1:1:nDof
            dNdx(a,i) = 0.;
            for j = 1:1:nDof
                dNdx(a,i) = dNdx(a,i) + dNdxi(a,j)*dxidx(j,i);
            end
        end
    end
    
    %Compute the element conductivity
    kel = kel + k*(dNdx(:,1)*dNdx(:,1)'+dNdx(:,2)*dNdx(:,2)')*w(intpt)*dt;
%     for a = 1:nNoEl
%         for b = 1:nNoEl
%             row = a;
%             col = b;
%             kel(row,col) = kel(row,col) + k*(dNdx(b,1)*dNdx(a,1)+dNdx(b,2)*dNdx(a,2))*w(intpt)*dt;
%         end
%     end
end
end