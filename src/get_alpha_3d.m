function alpha = get_alpha(tri,pb,T1p,Phi)

% compute alpha matrix by collocation

t = tri.ConnectivityList;
n = max(pb(:,2));
alpha = sparse(n,n);
Te = pointLocation(tri,T1p);        % numbers of elements the points T1p are located in
%I = find(isnan(Te)); size(I)
%scatter(p(I,1),p(I,2),20,'filled');

for k = 1:n                         % loop over all points
    e = Te(k);
    Tns = t(e,1:4);                 % node numbers of element Te(k)
    alpha(k,pb(Tns,2)') = (Phi{e}*[1 T1p(k,:)]')';
end
