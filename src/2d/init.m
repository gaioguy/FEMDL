% various little functions

addpath('../../src')
clear all

% evaluation at end time
at_tf = @(A) squeeze(A(:,:,end,:));

% symmetric part of matrix
sym = @(A) 0.5*(A+A');

% derivative of f at x in direction d
df = @(f,x,d) imag(f(x + i*eps*d).')/eps; 

% Jacobi matrices of f at the rows of x
D = @(f,x) permute(cat(3,df(f,x,[1,0]),df(f,x,[0,1])),[1 3 2]); 

mr = @(x) reshape(x,2*size(x,1),1);
rm = @(x) reshape(x,size(x,1)/2,2);
