function [c,minipr,minlevel,iprs] = curve_of_optimal_ipr(mesh,T,u,lb,ub)

n = size(u,1);
p = mesh.p;
p = p(1:n,:);
tri = alphaShape(p(:,1),p(:,2),1);
t = alphaTriangulation(tri);
lv = linspace(lb,ub,200);
iprs = ones(length(lv),1)*inf;
for k = 1:length(lv)
    ct = tricontourc(p,t,u,[lv(k) lv(k)])';
    if ~isempty(ct)
        c = [ct(1).X' ct(1).Y'];  
%        Tc = T(c);
        Tc = mod(T(c)+[pi,0], 2*pi); % hack for standard map
        if is_closed(c) & is_closed(Tc)
            iprs(k) = 0.5*(ipr(c) + ipr(Tc));
        end
        % subplot(2,1,1); plot(c(:,1),c(:,2),'ko-','linewidth',1); axis([0 1 0 1]);
        % subplot(2,1,2); plot(Tc(:,1),Tc(:,2),'ko-','linewidth',1); 
        % axis([0 1 0 1]); drawnow;
        % pause
    end
end
[minipr,K] = min(iprs);                       % index of optimal curve
minlevel = lv(K);
ct = tricontourc(p,t,u,[lv(K) lv(K)])';
c = [ct(1).X' ct(1).Y'];

