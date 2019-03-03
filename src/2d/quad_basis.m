function a = quad_basis(p,t,pb,h)

%% QUAD_BASIS compute integrals of basis functions
%
%  a = QUAD_BASIS(p,t,pb,h) computes the integrals of the P1 Lagrange basis
%  functions over the domain of the triangulation w.r.t. the density h
%   p: (n x 2), one node per row
%   t: (m x 3), integers, each row defines a triangle by indexing into p
%   pb: (n x 2), node pb(i,2) maps to pb(i,1) (for perodic boundaries)
%   h: function handle, density for integral
%
% based on code from ifem by Long Chen
%
% (C) 2017 by O. Junge, see COPYRIGHT 

n = max(pb(:,2));                             % number of nodes
m = size(t,1);   
deg = 3; 

A = zeros(m,3); 
[lam,w] = quadpts(deg);
[~,area] = gradbasis(p,t);                    % compute areas of triangles
phi = lam;                                    % basis functions at quadrature points
nq = size(lam,1);
for k = 1:nq  
    phi_p = [phi(k,1) phi(k,2) phi(k,3)];     % values of basis functions at quadrature point
    x = lam(k,1)*p(t(:,1),:) ...              % quadrature points in the x-y coordinate
      + lam(k,2)*p(t(:,2),:) ...
      + lam(k,3)*p(t(:,3),:);
    A = A + w(k)*h(x)*phi_p;
end
A = A.*area;
a = accumarray(pb(t(:),2),A(:),[n 1]);        % a(i) = \int \phi_i h(x) dx

