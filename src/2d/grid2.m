function nodes = grid(nx, ny)

x1 = linspace(0,1,nx); 
y1 = linspace(0,1,ny);
[X,Y] = meshgrid(x1,y1); 
nodes = [X(:) Y(:)];
