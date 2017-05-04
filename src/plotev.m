function plotev2(t,p,pb,v,mesh)

w(pb(:,1)) = v(pb(:,2)); 
w = w/norm(w,inf); 
trisurf(t,p(:,1),p(:,2),zeros(size(p,1),1),w); 
shading interp;
if mesh,
    hold on; triplot(t,p(:,1),p(:,2),'k','linewidth',0.1);
end
view(2), axis equal, axis tight, caxis([-1,1]),
