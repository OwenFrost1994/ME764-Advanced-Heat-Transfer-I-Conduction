function [g]=gmake(M,N,dx,dy,gdd)
g = zeros(M*N,1);
%internal nodes
for i=2:1:M-1
    for j=2:1:N-1
        g(M*(j-1)+i) = gdd*dx*dy;
    end
end
%edge nodes
for i=2:1:M-1
    g(M*(1-1)+i) = gdd*dx*dy/2;
    g(M*(N-1)+i) = gdd*dx*dy/2;
end
for j=2:1:N-1
    g(M*(j-1)+1) = gdd*dx*dy/2;
    g(M*(j-1)+M) = gdd*dx*dy/2;
end
%corner nodes
g(1) = gdd*dx*dy/4;
g(M) = gdd*dx*dy/4;
g(M*(N-1)+1) = gdd*dx*dy/4;
g(M*N) = gdd*dx*dy/4;
end