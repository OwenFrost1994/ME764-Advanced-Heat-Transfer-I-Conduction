function [nNoEl]=numNoEl(nDof,order)
if(order==1)
    nNoEl=4;    %% number of Node in an 2D Element for order 1(Q4)
end
if(order==2)
    nNoEl=8;     %% number of Node in an 2D Element for order 2(Q8)
end
end