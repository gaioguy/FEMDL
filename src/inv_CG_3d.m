function CG = inv_CG_3d(v,x,tspan)

%% INV_CG_3D compute averaged inverse Cauchy-Green tensor
%
% CG = inv_CG_3d(v,x,tspan)
%   v: velocity field, function handle of the form v(t,x)
%   x: (n x 3), intial data in 3D
%   tspan: vector of time instances
% returns a (n x 6) matrix CG with CG(i,:) = [G(1,1) G(1,2) G(1,3) G(2,2) G(2,3) G(3,3)]
% where G(i,:) is the averaged inverse Cauch-Green tensor for initial condition
% x(i,:).
%
% inspired by code by Alireza Hadjighasem, alirezah@mit.edu
% with contributions by Daniel Karrasch, karrasch@ma.tum.de
% and Andreas Denner
%
% (C) 2017 by O. Junge and G. Froyland, see COPYRIGHT 

xi = x(:,1); yi = x(:,2); zi = x(:,3);
rho.x = 1e-8; rho.y = 1e-8; rho.z = 1e-8;

q = numel(tspan);
m = numel(xi);
Nrad = 6;
x0 = repmat(xi,1,Nrad);
y0 = repmat(yi,1,Nrad);
z0 = repmat(zi,1,Nrad);
%%
x0(:,1)=xi+rho.x;
x0(:,3)=xi-rho.x;
y0(:,2)=yi+rho.y;
y0(:,4)=yi-rho.y;
z0(:,5)=zi+rho.z;
z0(:,6)=zi-rho.z;

%% time integration
[xt,yt,zt] = integrator3(v,x0(:),y0(:),z0(:),tspan);
xt = reshape(xt,q, m, Nrad); % Number of time steps x number of pointsxdimension x 6
yt = reshape(yt,q, m, Nrad);
zt = reshape(zt,q, m, Nrad);

%% finite-differencing on a time-space grid
F11 = (xt(:,:,1)-xt(:,:,3))/(2*rho.x); 
F12 = (xt(:,:,2)-xt(:,:,4))/(2*rho.x); 
F13 = (xt(:,:,5)-xt(:,:,6))/(2*rho.x); 

F21 = (yt(:,:,1)-yt(:,:,3))/(2*rho.y);
F22 = (yt(:,:,2)-yt(:,:,4))/(2*rho.y);
F23 = (yt(:,:,5)-yt(:,:,6))/(2*rho.y);

F31 = (zt(:,:,1)-zt(:,:,3))/(2*rho.z);
F32 = (zt(:,:,2)-zt(:,:,4))/(2*rho.z);
F33 = (zt(:,:,5)-zt(:,:,6))/(2*rho.z);

%% Compute inverse of Cauchy-Green tensor

F = zeros(3,3,numel(F11));
F(1,1,:) = F11(:); F(1,2,:) = F12(:); F(1,3,:) = F13(:);
F(2,1,:) = F21(:); F(2,2,:) = F22(:); F(2,3,:) = F23(:);
F(3,1,:) = F31(:); F(3,2,:) = F32(:); F(3,3,:) = F33(:);

IF = permute(MultiSolver(F,eye(3)),[1 3 2]);
IF11 = reshape(IF(1,1,:),q,m);
IF12 = reshape(IF(1,2,:),q,m);
IF13 = reshape(IF(1,3,:),q,m);
IF21 = reshape(IF(2,1,:),q,m);
IF22 = reshape(IF(2,2,:),q,m);
IF23 = reshape(IF(2,3,:),q,m);
IF31 = reshape(IF(3,1,:),q,m);
IF32 = reshape(IF(3,2,:),q,m);
IF33 = reshape(IF(3,3,:),q,m);

%% CG'*CG
IC11 = IF11.^2+IF12.^2+IF13.^2;
IC22 = IF21.^2+IF22.^2+IF23.^2;
IC33 = IF31.^2+IF32.^2+IF33.^2;
IC12 = IF11.*IF21+IF12.*IF22+IF13.*IF23;
IC13 = IF11.*IF31+IF12.*IF32+IF13.*IF33;
IC23 = IF21.*IF31+IF22.*IF32+IF23.*IF33;

%% average over time steps
meanIC11 = 1/q*sum(IC11,1);
meanIC22 = 1/q*sum(IC22,1);
meanIC33 = 1/q*sum(IC33,1);
meanIC12 = 1/q*sum(IC12,1);
meanIC13 = 1/q*sum(IC13,1);
meanIC23 = 1/q*sum(IC23,1);

%% output
CG = [meanIC11(:) meanIC12(:) meanIC13(:) meanIC22(:) meanIC23(:) meanIC33(:)];