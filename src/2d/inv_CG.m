function [CG,DT] = inv_CG(v,x,tspan)

%% INV_CG compute averaged Cauchy-Green tensor
%
% CG = inv_CG(v,x,tspan)
%   v: velocity field, function handle of the form v(t,x)
%   x: (n x 2), intial data in 2D
%   tspan: vector of time instances
% returns a (n x 3) matrix CG with CG(i,:) = [G(1,1) G(1,2) G(2,2)]
% and G the averaged inverse Cauch-Green tensor for initial condition
% x(i,:)
%
% inspired by code by Alireza Hadjighasem, alirezah@mit.edu
% with contributions by Daniel Karrasch, karrasch@ma.tum.de
%
% (C) 2017 by O. Junge and G. Froyland, see COPYRIGHT 

xi = x(:,1); yi = x(:,2);
rho.x = 1e-8; rho.y = rho.x;

%% now for the Jacobian etc.
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
F11 = (xt(:,:,1)-xt(:,:,3))/(2*rho.x);
F12 = (xt(:,:,2)-xt(:,:,4))/(2*rho.x);
F21 = (yt(:,:,1)-yt(:,:,3))/(2*rho.y);
F22 = (yt(:,:,2)-yt(:,:,4))/(2*rho.y);

%% getting the components of inv(DT)
IF11 = arrayfun(@(a,b,c,d) [1, 0]*([a, b; c, d]\[1; 0]),F11,F12,F21,F22);
IF12 = arrayfun(@(a,b,c,d) [1, 0]*([a, b; c, d]\[0; 1]),F11,F12,F21,F22);
IF21 = arrayfun(@(a,b,c,d) [0, 1]*([a, b; c, d]\[1; 0]),F11,F12,F21,F22);
IF22 = arrayfun(@(a,b,c,d) [0, 1]*([a, b; c, d]\[0; 1]),F11,F12,F21,F22);

IC1 = IF11.^2+IF12.^2;
IC2 = IF11.*IF21+IF12.*IF22;
IC3 = IF21.^2+IF22.^2;

meanIC1 = 1/q*sum(IC1,1);
meanIC2 = 1/q*sum(IC2,1);
meanIC3 = 1/q*sum(IC3,1);

CG = [meanIC1(:) meanIC2(:) meanIC3(:)];

meanDT1 = 1/q*sum(F11,1);
meanDT2 = 1/q*sum(F12,1);
meanDT3 = 1/q*sum(F21,1);
meanDT4 = 1/q*sum(F22,1);

DT = [meanDT1(:) meanDT2(:) meanDT3(:) meanDT4(:)];

