function Q = dotA3_sym(x,A,y)

% general scalar product in R^3

% A is m x 6, each row represents m symmetric 3 x 3 matrices with
% indices [1 2 3; 
%          2 4 5;
%          3 5 6]
% x,y are m x 3, each row is a 3-vector x(m), y(m)
% Q is m x 1, each row is the dot product x(m)'*A(m)*y(m)

Q = x(:,1).*(A(:,1).*y(:,1) + A(:,2).*y(:,2) + A(:,3).*y(:,3)) ...
  + x(:,2).*(A(:,2).*y(:,1) + A(:,4).*y(:,2) + A(:,5).*y(:,3)) ...
  + x(:,3).*(A(:,3).*y(:,1) + A(:,5).*y(:,2) + A(:,6).*y(:,3));