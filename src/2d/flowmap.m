function y = flowmap(v,x,tspan,tol)

nt = numel(tspan);
[nx, dim] = size(x);
y = zeros(nx,dim,nt);

if nargin < 4
    tol = 1e-3;
end
options = odeset('RelTol',tol,'AbsTol',tol);
[~,F] = ode45(v,tspan,[x(:,1); x(:,2)],options);
[nF,mF] = size(F);

if nt==2, 
    ts = [1 nF]; 
else
    ts = 1:nt;
end

for k = 1:length(ts)
    y(:,1,k) = F(ts(k),1:mF/2);
    y(:,2,k) = F(ts(k),mF/2+1:mF);
end
