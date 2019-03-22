function m = trimesh(nx, ny, type)

n = nx*ny;

if type=='grid'
    x1 = linspace(xmin,xmax,nx); 
    y1 = linspace(ymin,ymax,ny);
    [X,Y] = meshgrid(x1,y1); 
    p = [X(:) Y(:)];
else
    error('implement me.');
end

if strcmp(top,'square')
    pb = [1:n; 1:n]';
elseif strcmp(top,'cylinder')
    pb = [1:n; [1:((nx-1)*ny), 1:ny]]'; 
else
    error('implement me.');
end    

Tri = alphaShape(p(:,1),p(:,2)); 
t = alphaTriangulation(Tri);
b = unique(boundaryFacets(Tri));

m.p = p;
m.t = t;
m.b = b;
m.pb = pb;
