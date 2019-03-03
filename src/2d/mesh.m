function mesh = mesh(ll, ur, ns, type, top)

xmin = ll(1); xmax = ur(1);
ymin = ll(2); ymax = ur(2);
nx = ns(1); ny = ns(2); n = nx*ny;

if type=='uniform'
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

mesh.p = p;
mesh.t = t;
mesh.b = b;
mesh.pb = pb;
