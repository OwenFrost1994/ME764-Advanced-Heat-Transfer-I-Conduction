%====================== No. integration points =============================
%
%   Defines the number of integration points:be used for
%   each element type
%
function n = numberofintegrationpoints(nDof,nNoEl)
if (nDof == 1) 
     n = nNoEl;
elseif (nDof == 2)
    if (nNoEl == 4)
        n = 4;
    end
    if (nNoEl == 8)
        n = 9;
    end
elseif (nDof == 3)
    if (nNoEl == 4)
        n = 1 ;
    end
    if (nNoEl == 8)
        n = 8;
    end
    if (nNoEl == 10)
        n = 4;
    end
    if (nNoEl == 20)
        n = 27;
    end
end
end