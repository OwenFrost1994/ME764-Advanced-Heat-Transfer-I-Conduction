function [Hin] = Hinmake(M,N,hinl,hinv,dx,dy)
if mod(N,2) == 0
    Hin = zeros(M*N,M*N);
    %internal nodes
    for i=2:1:M-1
        for j=2:1:N/2
            Hin(M*(j-1)+i,M*(j-1)+i) = hinl*dx*dy;
        end
        for j=N/2+1:1:N-1
            Hin(M*(j-1)+i,M*(j-1)+i) = hinv*dx*dy;
        end
    end
    %edge nodes
    for i=2:1:M-1
        Hin(M*(1-1)+i,M*(1-1)+i) = hinl*dx*dy/2;
        Hin(M*(N-1)+i,M*(N-1)+i) = hinv*dx*dy/2;
    end
    for j=2:1:N/2
        Hin(M*(j-1)+1,M*(j-1)+1) = hinl*dx*dy/2;
        Hin(M*(j-1)+M,M*(j-1)+M) = hinl*dx*dy/2;
    end
    for j=N/2+1:1:N-1
        Hin(M*(j-1)+1,M*(j-1)+1) = hinv*dx*dy/2;
        Hin(M*(j-1)+M,M*(j-1)+M) = hinv*dx*dy/2;
    end
    %corner nodes
    Hin(1,1) = hinl*dx*dy/4;
    Hin(M,M) = hinl*dx*dy/4;
    Hin(M*(N-1)+1,M*(N-1)+1) = hinv*dx*dy/4;
    Hin(M*N,M*N) = hinv*dx*dy/4;
else
    Hin = zeros(M*N,M*N);
    %internal nodes
    for i=2:1:M-1
        for j=2:1:floor(N/2)
            Hin(M*(j-1)+i,M*(j-1)+i) = hinl*dx*dy;
        end
        for j=floor(N/2)+2:1:N-1
            Hin(M*(j-1)+i,M*(j-1)+i) = hinv*dx*dy;
        end
    end
    %edge nodes
    for i=2:1:M-1
        Hin(M*(1-1)+i,M*(1-1)+i) = hinl*dx*dy/2;
        Hin(M*(N-1)+i,M*(N-1)+i) = hinv*dx*dy/2;
        Hin(M*(floor(N/2))+i,M*(floor(N/2))+i) = (hinl*dx*dy+hinv*dx*dy)/2;
    end
    for j=2:1:floor(N/2)
        Hin(M*(j-1)+1,M*(j-1)+1) = hinl*dx*dy/2;
        Hin(M*(j-1)+M,M*(j-1)+M) = hinl*dx*dy/2;
    end
    for j=floor(N/2)+2:1:N-1
        Hin(M*(j-1)+1,M*(j-1)+1) = hinv*dx*dy/2;
        Hin(M*(j-1)+M,M*(j-1)+M) = hinv*dx*dy/2;
    end
    %corner nodes
    Hin(1,1) = hinl*dx*dy/4;
    Hin(M,M) = hinl*dx*dy/4;
    Hin(M*(N-1)+1,M*(N-1)+1) = hinv*dx*dy/4;
    Hin(M*N,M*N) = hinv*dx*dy/4;
end
end