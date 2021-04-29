function plot_lcs(X,Y,Z,I,val,col,alpha)

[n,n,n] = size(X);
J = find(I~=val); length(J);
I(J) = 1000;
ww = reshape(I,n,n,n);
patch(isosurface(X,Y,Z,ww,val),'FaceColor',col,'EdgeColor','none','FaceAlpha',alpha);
patch(isocaps(X,Y,Z,ww,val,'below'),'FaceColor',col,'EdgeColor','none','FaceAlpha',alpha)


