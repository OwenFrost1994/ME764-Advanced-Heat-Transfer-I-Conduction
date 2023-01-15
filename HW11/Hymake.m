function [H]=Hymake(M,N,hfd,dhdev,beta,dx,dy,L)
H = zeros(M*N,M*N);
%edge nodes
for j=2:1:N-1
    H(M*(j-1)+1,M*(j-1)+1) = (hy(dy*(j-2),hfd,dhdev,beta,L)+hy(dy*(j),hfd,dhdev,beta,L))/2*dy;
end
%corner nodes
H(1,1) = (hy(0,hfd,dhdev,beta,L))*dy/2;
H(M*(N-1)+1,M*(N-1)+1) = (hy((N-1)*dy,hfd,dhdev,beta,L))*dy/2;
end