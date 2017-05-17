function dz = identity(t,z)

n = numel(z)/2;
dz(1:n,1) = zeros(n,1);
dz(n+1:2*n,1) = zeros(n,1);



