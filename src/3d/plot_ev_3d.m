function plot_ev_3d(X,Y,Z,w,theta)

[n,n,n] = size(X);
w = w/norm(w,inf); 
ww = reshape(w, n, n, n);

h1 = patch(isosurface(X,Y,Z,ww,theta),'FaceAlpha',0.5);
isonormals(X,Y,Z,ww,h1); 
h1.FaceColor = 'red'; h1.EdgeColor = 'none';
patch(isocaps(X,Y,Z,ww,theta,'above'),'FaceColor','interp','EdgeColor','none','FaceAlpha',0.5)

h2 = patch(isosurface(X,Y,Z,ww,-theta),'FaceAlpha',0.5);
isonormals(X,Y,Z,ww,h2);
h2.FaceColor = 'blue'; h2.EdgeColor = 'none'; 
patch(isocaps(X,Y,Z,ww,-theta,'below'),'FaceColor','interp','EdgeColor','none','FaceAlpha',0.5)
patch(isocaps(X,Y,Z,ww,theta,'below'),'FaceColor','interp','EdgeColor','none','FaceAlpha',0.25)
