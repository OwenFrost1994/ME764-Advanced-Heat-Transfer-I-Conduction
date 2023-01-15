function [C]=Ctmake(M,N,rho,dx,dy,t)
C = zeros(M*N,M*N);
%internal nodes
for i=2:1:M-1
    for j=2:1:N-1
        C(M*(j-1)+i,M*(j-1)+i) = rho*ct(t(M*(j-1)+i))*dx*dy;
    end
end
%edge nodes
for i=2:1:M-1
    C(M*(1-1)+i,M*(1-1)+i) = rho*ct(t(M*(1-1)+i))*dx*dy/2;
    C(M*(N-1)+i,M*(N-1)+i) = rho*ct(t(M*(N-1)+i))*dx*dy/2;
end
for j=2:1:N-1
    C(M*(j-1)+1,M*(j-1)+1) = rho*ct(t(M*(j-1)+1))*dx*dy/2;
    C(M*(j-1)+M,M*(j-1)+M) = rho*ct(t(M*(j-1)+M))*dx*dy/2;
end
%corner nodes
C(1,1) = rho*ct(t(1))*dx*dy/4;
C(M,M) = rho*ct(t(M))*dx*dy/4;
C(M*(N-1)+1,M*(N-1)+1) = rho*ct(t(M*(N-1)+1))*dx*dy/4;
C(M*N,M*N) = rho*ct(t(M*N))*dx*dy/4;
end