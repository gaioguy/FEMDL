function plotf(m,f,plot_edges)

%% PLOTF plot scalar function on triangulation
%
% PLOTF(t,p,pb,F,mesh) plots the scalar function given by 
%   m: triangle mesh as produced by trimesh
%   f: (n x 1), values of the f on the nodes of m
%   edges: logical, determines whether triangulation is plotted as well
%
% (C) 2019 by O. Junge, see COPYRIGHT 

p = m.p; t = m.t; pb = m.pb;
w(pb(:,1)) = f(pb(:,2)); 
trisurf(t,p(:,1),p(:,2),zeros(size(p,1),1),w);  
shading interp;
if plot_edges,
   hold on; triplot(t,p(:,1),p(:,2),'k','linewidth',0.1);
end
view(2), axis equal, axis tight
mf = max(abs(f)); caxis([-mf mf]);
drawnow
