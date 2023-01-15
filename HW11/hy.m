function [h]=hy(y,hfd,dhdev,beta,L)
h = hfd+dhdev/(1+beta*y/L);
end