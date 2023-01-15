function [tinfin, tinfout] = tinfmake(M,N,tin,tout)
tinfin = tin*ones(M*N,1);
tinfout = tout*ones(M*N,1);
end