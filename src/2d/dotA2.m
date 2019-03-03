function Q = dotA2(x,A,y)

% general scalar product in R^2

% A is m x 4, each row represents a 2 x 2 matrix A(m) = [a_11 a_12 a_21 a_22]
% x,y are m x 2, each row is a 2-vector x(m), y(m)
% Q is m x 1, each row is the dot product x(m)'*A(m)*y(m)

Q = x(:,1).*(A(:,1).*y(:,1) + A(:,2).*y(:,2)) ...
  + x(:,2).*(A(:,3).*y(:,1) + A(:,4).*y(:,2));