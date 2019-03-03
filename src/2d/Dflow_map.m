function DT = Dflow_map(v,x,tspan)

%% DFLOW_MAP compute Jacobi matrices of flow map by finite differencing
%
% DT = Dflow_map(v,x,tspan)
%   v: velocity field, function handle of the form v(t,x)
%   x: (n x 2), intial data in 2D
%   tspan: q-vector of time instances
% returns a (q x m x 2 x 2) tensor DT with DT(t,j,:,:) the 2 x 2 Jacobi
% matrix at time t for initial condition x(i,:)
%
% inspired by code by Alireza Hadjighasem, alirezah@mit.edu
%
% (C) 2018 by O. Junge, see COPYRIGHT 

xi = x(:,1); yi = x(:,2);
rho.x = 1e-8; rho.y = rho.x;

q = numel(tspan);
m = numel(xi);
Nrad = 4;
xt = zeros(m,Nrad);
yt = zeros(m,Nrad);
for k = 1:Nrad
    xt(:,k) = xi + rho.x*cos( (k-1)*pi/2 );
    yt(:,k) = yi + rho.y*sin( (k-1)*pi/2 );
end

%% Time integration
[xt,yt] = integrator(v,xt(:),yt(:),tspan);

xt = reshape(xt,q, m, Nrad);
yt = reshape(yt,q, m, Nrad);

%% finite-differencing on a time-space grid
DT = zeros(2,2,q,m);
DT(1,1,:,:) = (xt(:,:,1)-xt(:,:,3))/(2*rho.x);
DT(1,2,:,:) = (xt(:,:,2)-xt(:,:,4))/(2*rho.x);
DT(2,1,:,:) = (yt(:,:,1)-yt(:,:,3))/(2*rho.y);
DT(2,2,:,:) = (yt(:,:,2)-yt(:,:,4))/(2*rho.y);