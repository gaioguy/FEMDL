function V = eval_p1(X,Y,Z,v,x)

%% EVAL_P1_3D evaluate piecewise linear functions
%
% V = EVAL_P1_3D(p,v,x) evaluates the piecewise linear functions defined
% on the nodes p with values V(:,i)
%   p: (n x 3), one node per row
%   v: (n x j), j function values on each node 
%   x: (m x 3). points, on which the functions are to be evaluated
% V(:,j) are the values of v(:,j) in the points x
%
% (C) 2017 by O. Junge, see COPYRIGHT 

k = size(v,2);
P = [2 1 3];
X = permute(X, P);
Y = permute(Y, P);
Z = permute(Z, P);
for j = 1:k,
    w = reshape(v(:,j),size(X));
    w = permute(w, P);
    vI{j} = griddedInterpolant(X,Y,Z,w);
    V(:,j) = vI{j}(x(:,1),x(:,2),x(:,3)); 
end
