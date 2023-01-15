function [H] = Hmake(M,N,h,dx,dy)
H = zeros(M*N,M*N);
%edge nodes
for i=2:1:M-1
    H(M*(N-1)+i,M*(N-1)+i) = h*dx;
end
for j=2:1:N-1
    H(M*(j-1)+1,M*(j-1)+1) = h*dy;
end
%corner nodes
H(1,1) = h*dy/2;
H(M*(N-1)+1,M*(N-1)+1) = h*(dx+dy)/2;
H(M*N,M*N) = h*dx/2;
end