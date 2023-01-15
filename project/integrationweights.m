%
%================= INTEGRATION WEIGHTS ==================================
%
%   Defines integration weights w_i
%
function w=integrationweights(nDof,nNoEl,npoints)

w=zeros(npoints,1);
%
%  1D elements
%
if (nDof == 1) 
     if (npoints == 1)
       w(1) = 2.;
     elseif (npoints == 2) 
       w = [1.,1.];
     elseif (npoints == 3) 
       w = [0.555555555,0.888888888,0.555555555];
     end
end
%
%  2D elements
%
if (nDof == 2) 
%
%    Rectangular element
%                  
if ( nNoEl==4 || nNoEl==8 ) 
    if (npoints == 1) 
        w(1) = 4.;
    elseif (npoints == 4) 
            w = [1.,1.,1.,1.];
    elseif (npoints == 9 ) 
        w1D = [0.555555555,0.888888888,0.55555555555];
        for j = 1:3
            for i = 1:3
                n = 3*(j-1)+i;
                w(n) = w1D(i)*w1D(j);
            end
        end
    end
end
end
    
if (nDof == 3)
%
%  3D elements
%
    if ( nNoEl==8 || nNoEl==20 )
        if (npoints == 1) 
            w(1) = 8.;
        elseif (npoints == 8) 
            w = [1.,1.,1.,1.,1.,1.,1.,1.];
        elseif (npoints == 27) 
            w1D = [0.555555555,0.888888888,0.55555555555];
            for k = 1:3
                for j = 1:3
                    for i = 1:3
                        n = 9*(k-1)+3*(j-1)+i;
                        w(n) = w1D(i)*w1D(j)*w1D(k);
                    end
                end
            end
        end
    end
end
end
