function [th]=that(s,x,beta)
th = exp(-x)./(s.*(s-1))-(beta+1).*exp(-sqrt(s).*x)./(s.*(s-1).*(beta+sqrt(s)));
end