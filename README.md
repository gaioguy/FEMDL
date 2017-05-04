# FEMDL

A MATLAB package for a finite element based discretization of dynamic Laplacians

## Usage
Define the vector field, the times at which you wish to evaluate the flow map and then the flow map itself.
```Matlab
vf = @(t,x) rotating_double_gyre_vf(t,x);
t0 = 0; tf = 1; nt = 2;  ts = linspace(t0,tf,nt);
T = @(x) flow_map(v,x,ts);
```
Next generate a set of nodes and integrate them by the flow map. The matrix pb serves to define periodic boundaries (here, the boundary is non-periodic).
```Matlab
n = 625; p0 = rand(n,2);                            
P = T(p0);  
for k = 1:nt, 
   p{k}=P(:,[k k+nt]);
end; 
pb = [1:n; 1:n]';
```
Then you are ready to construct the stiffness and mass matrix:
```Matlab
D = sparse(n,n); 
for k = 1:nt                                        
    t{k} = delaunay(p{k});                          
    [Dt,Mt] = assemble(p{k},t{k},pb,kron([1 0 1],ones(size(t{k},1),1)));
    D = D + Dt/nt;  M = M + Mt/nt; 
end;
```
Next solve the corresponding eigenproblem
```Matlab
[V,L] = eigs(D,M,10,'SM');                          
[lam,I] = sort(diag(L),'descend'); V = V(:,I);
```
and use the eigenvectors to construct a partition of coherent sets by kmeans clustering:
```Matlab
n1 = 200; x1 = linspace(0,1,n1); [X1,Y1] = meshgrid(x1,x1); 
W = eval_p1(p{1},V(:,1:3),[X1(:) Y1(:)]);     
idx = kmeans(W, size(V1,2));               
```
Finally, plot the spectrum
```Matlab
plot(lam,'s','markerfacecolor','b'); 
```
some eigenvector
```Matlab
plotev(t{1},p{1},pb,V(:,2),1); colorbar
```
and the partition
```Matlab
surf(X1,Y1,reshape(idx,n1,n1)); shading flat
view(2); axis equal; axis tight; colorbar
colormap(70*[ 1 1 1; 2 2 2; 3 3 3]/255);
```

There are more demo computations provided in experiments/ which should be
executed after changing into the corresponding folder.

## Reference
G. Froyland and O. Junge. *Robust FEM-based extraction of finite-time 
coherent sets using scattered, sparse, and incomplete trajectories*.
arXiv, 2017.


