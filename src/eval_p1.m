function V = eval_p1(p,v,x)

%% EVAL_P1 evaluate piecewise linear functions
%
% V = EVAL_P1(p,v,x) evaluates the piecewise linear functions defined
% on the nodes p with values V(:,i)
%   p: (n x 2), one node per row
%   v: (n x j), j function values on each node 
%   x: (m x 2). points, on which the functions are to be evaluated
% V(:,j) are the values of v(:,j) in the points x
%
% (C) 2017 by O. Junge and G. Froyland, see COPYRIGHT 

k = size(v,2);
for j = 1:k,
    vI{j} = scatteredInterpolant(p(:,1),p(:,2),v(:,j));
    V(:,j) = vI{j}(x(:,1),x(:,2)); 
end
