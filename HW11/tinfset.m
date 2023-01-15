function [Tinf] = tinfset(M,N,tinf)
Tinf = zeros(M*N,1);
%edge nodes
for j=2:1:N-1
    Tinf(M*(j-1)+1) = tinf;
end
%corner nodes
Tinf(1) = tinf;
Tinf(M*(N-1)+1) = tinf;
end