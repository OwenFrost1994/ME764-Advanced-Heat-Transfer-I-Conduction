function [Tinf] = tinfmake(M,N,tinf)
Tinf = zeros(M*N,1);
%edge nodes
for i=2:1:M-1
    Tinf(M*(N-1)+i) = tinf;
end
for j=2:1:N-1
    Tinf(M*(j-1)+1) = tinf;
end
%corner nodes
Tinf(1) = tinf;
Tinf(M*(N-1)+1) = tinf;
Tinf(M*N) = tinf;
end