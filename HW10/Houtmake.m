function [Hout]=Houtmake(M,N,hout,dx,dy)
Hout = zeros(M*N,M*N);
%internal nodes
for i=2:1:M-1
    for j=2:1:N-1
        Hout(M*(j-1)+i,M*(j-1)+i) = hout*dx*dy;
    end
end
%edge nodes
for i=2:1:M-1
    Hout(M*(1-1)+i,M*(1-1)+i) = hout*dx*dy/2;
    Hout(M*(N-1)+i,M*(N-1)+i) = hout*dx*dy/2;
end
for j=2:1:N-1
    Hout(M*(j-1)+1,M*(j-1)+1) = hout*dx*dy/2;
    Hout(M*(j-1)+M,M*(j-1)+M) = hout*dx*dy/2;
end
%corner nodes
Hout(1,1) = hout*dx*dy/4;
Hout(M,M) = hout*dx*dy/4;
Hout(M*(N-1)+1,M*(N-1)+1) = hout*dx*dy/4;
Hout(M*N,M*N) = hout*dx*dy/4;
end