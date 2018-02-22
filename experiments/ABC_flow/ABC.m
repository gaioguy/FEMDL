function dX = ABC(t,X)

n = numel(X);
x = X(1:n/3);
y = X(n/3+1:2/3*n);
z = X(2/3*n+1:end);

A = sqrt(3); B = sqrt(2); C = 1;
A = A*ones(n/3,1);
A = A + 0.5*t.*sin(pi*t);

dx = A.*sin(2*pi*z)+C.*cos(2*pi*y);
dy = B.*sin(2*pi*x)+A.*cos(2*pi*z);
dz = C.*sin(2*pi*y)+B.*cos(2*pi*x);

dX = [dx;dy;dz];