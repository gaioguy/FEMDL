function plotf(p,t,pb,f,mesh)

%% PLOTF plot scalar function on triangulation
%
% PLOTF(t,p,pb,v,mesh) plots the scalar function given by 
%   p: (n x 2), one node per row
%   t: (m x 3), integers, each row defines a triangle by indexing into p
%   pb: (n x 2), node pb(i,2) maps to pb(i,1) (for perodic boundaries)
%   f: (n x 1), values of the function on the nodes
%   mesh: logical, determines whether triangulation is plotted as well
%
% (C) 2017 by O. Junge and G. Froyland, see COPYRIGHT 

w(pb(:,1)) = f(pb(:,2)); 
trisurf(t,p(:,1),p(:,2),zeros(size(p,1),1),w); 
shading interp;
if mesh,
    hold on; triplot(t,p(:,1),p(:,2),'k','linewidth',0.1);
end
view(2), axis equal, axis tight
