%====================== INTEGRATION POINTS ==================================
%
%   Defines positions of integration points
%
function xi=integrationpoints(nDof,nNoEl,npoints)

xi=zeros(nDof,npoints);
%
%  1D elements
%
if (nDof == 1) 
    if (npoints==1) 
        xi(1,1) = 0.;
    elseif (npoints == 2) 
        xi(1,1) = -0.5773502692;
        xi(1,2) = -xi(1,1);
    elseif (npoints == 3) 
        xi(1,1) = -0.7745966692;
        xi(1,2) = 0.0;
        xi(1,3) = -xi(1,1);
    end
end
%
% 2D elements
%
if (nDof == 2)
%    Rectangular element
    if ( nNoEl==4 || nNoEl==8 ) 
        if (npoints == 1) 
            xi(1,1) = 0.;
            xi(2,1) = 0.;
        elseif (npoints == 4) 
            xi(1,1) = -0.5773502692;
            xi(2,1) = xi(1,1);
            xi(1,2) = -xi(1,1);
            xi(2,2) = xi(1,1);
            xi(1,3) = xi(1,1);
            xi(2,3) = -xi(1,1);
            xi(1,4) = -xi(1,1);
            xi(2,4) = -xi(1,1);
        elseif (npoints == 9) 
            xi(1,1) = -0.7745966692;
            xi(2,1) = xi(1,1);
            xi(1,2) = 0.0;
            xi(2,2) = xi(1,1);
            xi(1,3) = -xi(1,1);
            xi(2,3) = xi(1,1);
            xi(1,4) = xi(1,1);
            xi(2,4) = 0.0;
            xi(1,5) = 0.0;
            xi(2,5) = 0.0;
            xi(1,6) = -xi(1,1);
            xi(2,6) = 0.0;
            xi(1,7) = xi(1,1);
            xi(2,7) = -xi(1,1);
            xi(1,8) = 0.;
            xi(2,8) = -xi(1,1);
            xi(1,9) = -xi(1,1);
            xi(2,9) = -xi(1,1);
        end
    end
%
%   3D elements
if (nDof == 3) 
%
%  3D elements
    if ( nNoEl==10 ) 
        if (npoints == 1)
            xi(1,1) = 0.25;
            xi(2,1) = 0.25;
            xi(3,1) = 0.25;
        elseif (npoints == 4)
            xi(1,1) = 0.58541020;
            xi(2,1) = 0.13819660;
            xi(3,1) = xi(2,1);
            xi(1,2) = xi(2,1);
            xi(2,2) = xi(1,1);
            xi(3,2) = xi(2,1);
            xi(1,3) = xi(2,1);
            xi(2,3) = xi(2,1);
            xi(3,3) = xi(1,1);
            xi(1,4) = xi(2,1);
            xi(2,4) = xi(2,1);
            xi(3,4) = xi(2,1);
        end
    elseif ( nNoEl==8 || nNoEl==20 )
        if (npoints == 1)
            xi(1,1) = 0.;
            xi(2,1) = 0.;
            xi(3,1) = 0.;
        elseif (npoints == 8)
            x1D = [-0.5773502692,0.5773502692];
            for k = 1:2
                for j = 1:2
                    for i = 1:2
                        n = 4*(k-1) + 2*(j-1) + i;
                        xi(1,n) = x1D(i);
                        xi(2,n) = x1D(j);
                        xi(3,n) = x1D(k);
                    end
                end
            end
        elseif (npoints == 27)
            x1D = [-0.7745966692,0.,0.7745966692];
            for k = 1:3
                for j = 1:3
                    for i = 1:3
                        n = 9*(k-1) + 3*(j-1) + i;
                        xi(1,n) = x1D(i);
                        xi(2,n) = x1D(j);
                        xi(3,n) = x1D(k);
                    end
                end
            end
        end
    end
end
end
