function plotf2(mesh, F, varargin)

%% PLOTF2 plot scalar functions on 2d triangulation
%
% PLOTF2(mesh,F) plots the scalar functions given by 
%   mesh:   mesh data structure
%   F: (n x k), values of the functions on the nodes
%   edges: logical, determines whether triangulation is plotted as well
%
% (C) 2018 by O. Junge, see COPYRIGHT 

ip = inputParser;
defaultLayout = [];
defaultEdges = 0;
defaultContour = 0;

addRequired(ip,'mesh');
addRequired(ip,'F');
addParamValue(ip,'layout',defaultLayout,@isvector);
addParamValue(ip,'edges',defaultEdges,@islogical);
addParamValue(ip,'contour',defaultContour,@islogical);

parse(ip, mesh, F, varargin{:});
layout  = ip.Results.layout;
edges   = ip.Results.edges;
contour = ip.Results.contour;

[n, m] = size(F);
p = mesh.p; t = mesh.t; pb = mesh.pb;

for k = 1:m
    f = F(:,k);
    if ~isempty(layout)
        subplot(layout(1),layout(2),k);
    end
    w(pb(:,1)) = f(pb(:,2)); 
    if contour
        tricontour(p,t,w,20);
    else
        trisurf(t,p(:,1),p(:,2),zeros(size(p,1),1),w);  shading interp;
        if edges,
            hold on; triplot(t,p(:,1),p(:,2),'k','linewidth',0.1);
        end
    end
    view(2), axis equal, axis tight
    mf = max(abs(f)); caxis([-mf mf]);
    colorbar
    drawnow
end
