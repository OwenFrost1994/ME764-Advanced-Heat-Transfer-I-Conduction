function [c]=ct(t)
if t>273.2
    c = 3500;
else if 260.2<t && t<273.2
        c = -41650000+312428*t-585.511*t^2;
    else
        c = 512.66+6.89*t;
    end
end
end