function [k]=kt(t)
if t>273.2
    k = 0.49;
else if 260.2<t && t<273.2
        k = 2.21-0.1331*(t-260.2);
    else
        k = 2135/(t^(1.235));
    end
end
end