% various little constants and functions

addpath('../../src')

getcolumns = @(A,I) A(:,I);
normed = @(v) v/norm(v,inf);

% sqrt of machine epsilon
eps2 = sqrt(eps);

% evaluation at end time
at_tf = @(A) squeeze(A(:,:,end,:));

% symmetric part of matrix
sym = @(A) 0.5*(A+A');

% derivative of f at x in direction d
df = @(f,x,d) imag(f(x + i*eps*d).')/eps; 
df_fd = @(f,x,d) (f(x+eps2*d)-f(x-eps2*d))'/abs(2*eps2);

% Jacobi matrices of f at the rows of x
D = @(f,x) permute(cat(3,df(f,x,[1,0]),df(f,x,[0,1])),[1 3 2]); 
D_fd = @(f,x) permute(cat(3,df_fd(f,x,[1,0]),df_fd(f,x,[0,1])),[1 3 2]); 

mr = @(x) reshape(x,2*size(x,1),1);
rm = @(x) reshape(x,size(x,1)/2,2);
