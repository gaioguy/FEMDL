function dx = ocean(t,x,U,V)

n = numel(x)/2;
dx = zeros(2*n,1);

dx(1:n,1)     = U( x(1:n,1),x(n+1:2*n,1),t*ones(n,1) );
dx(n+1:2*n,1) = V( x(1:n,1),x(n+1:2*n,1),t*ones(n,1) );
