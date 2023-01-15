function [K]=Kmake(M,N,k,dx,dy)
K = zeros(M*N,M*N);
%internal nodes
for i=2:1:M-1
    for j=2:1:N-1
        K(M*(j-1)+i,M*(j-1)+i-1) = -k*dy/dx;
        K(M*(j-1)+i,M*(j-1)+i+1) = -k*dy/dx;
        K(M*(j-1)+i,M*(j-2)+i) = -k*dx/dy;
        K(M*(j-1)+i,M*(j)+i) = -k*dx/dy;
        K(M*(j-1)+i,M*(j-1)+i) = -sum(K(M*(j-1)+i,:));
    end
end
%edge nodes
for i=2:1:M-1
    K(i,i-1) = -k*dy/dx/2;
    K(i,i+1) = -k*dy/dx/2;
    K(i,M+i) = -k*dx/dy;
    K(i,i) = -sum(K(i,:));
    K(M*(N-1)+i,M*(N-1)+i-1) = -k*dy/dx/2;
    K(M*(N-1)+i,M*(N-1)+i+1) = -k*dy/dx/2;
    K(M*(N-1)+i,M*(N-2)+i) = -k*dx/dy;
    K(M*(N-1)+i,M*(N-1)+i) = -sum(K(M*(N-1)+i,:));
end
for j=2:1:N-1
    K(M*(j-1)+1,M*(j)+1) = -k*dx/dy/2;
    K(M*(j-1)+1,M*(j-2)+1) = -k*dx/dy/2;
    K(M*(j-1)+1,M*(j-1)+2) = -k*dy/dx;
    K(M*(j-1)+1,M*(j-1)+1) = -sum(K(M*(j-1)+1,:));
    K(M*(j-1)+M,M*(j)+M) = -k*dx/dy/2;
    K(M*(j-1)+M,M*(j-2)+M) = -k*dx/dy/2;
    K(M*(j-1)+M,M*(j-1)+M-1) = -k*dy/dx;
    K(M*(j-1)+M,M*(j-1)+M) = -sum(K(M*(j-1)+M,:));
end
%corner nodes
K(1,2) = -k*dy/dx/2;
K(1,M+1) = -k*dx/dy/2;
K(1,1) = -sum(K(1,:));

K(M,M-1) = -k*dy/dx/2;
K(M,2*M) = -k*dx/dy/2;
K(M,M) = -sum(K(M,:));

K(M*(N-1)+1,M*(N-1)+2) = -k*dy/dx/2;
K(M*(N-1)+1,M*(N-2)+1) = -k*dx/dy/2;
K(M*(N-1)+1,M*(N-1)+1) = -sum(K(M*(N-1)+1,:));

K(M*N,M*N-1) = -k*dy/dx/2;
K(M*N,M*N-M) = -k*dx/dy/2;
K(M*N,M*N) = -sum(K(M*N,:));
end