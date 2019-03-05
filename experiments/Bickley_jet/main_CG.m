addpath('../../src','../../src/2d'); clear all; 
cols = @(A,I) A(:,I);
at_tf = @(A) squeeze(A(:,:,end,:));

%% flow map 
t0 = 0; days = 60*60*24; tf = 40*days; 
vf = @bickleyjet;
DT  = @(x) at_tf(Dflow_map(vf,x,[t0 tf]));  % space derivative of flow map
DL = @(DT) 0.5*(eye(2) + inv(DT)*inv(DT)'); % dynamic Laplace
DLx = @(x) fapply1(DL, DT(x));              % evaluate DL at each row of x

%% triangulation
nx = 100;  ny = floor(nx/20*6);  n = nx*ny;
[X,Y] = meshgrid(linspace(0,20,nx),linspace(-3,3,ny)); 
p = [X(:) Y(:)];
pb = [1:n; [1:((nx-1)*ny), 1:ny]]';         % boundary periodic in x
t = delaunay(p); 

%% assembly
deg = 1;                                    % degree of quadrature
tic; A = triquad(p,t,DLx,deg); toc          % integrate DL on triangles
tic; [D,M] = assemble2(p,t,pb,A); toc        % assemble stiffness and mass matrices

%% solve eigenproblem
tic; [V,L] = eigs(D,M,15,'SM'); toc
[lam,ord] = sort(diag(L),'descend'); V = V(:,ord);

%% plot spectrum
figure(1); clf; plot(lam,'*'); axis tight, axis square
xlabel('$k$'); ylabel('$\lambda_k$')

%% plot eigenvector
figure(2), plotf(p,t,pb,normed(V(:,3)),0); caxis([-1,1]); colorbar
ylabel('$y$');  xlabel('$x$'); 

%% compute partition
nx1 = 400; ny1 = nx1/20*6; x1 = linspace(0,20,nx1); y1 = linspace(-3,3,ny1);
[X1,Y1] = meshgrid(x1,y1); 
nc = 8;
V1 = eval_p1(p,V(pb(:,2),1:nc),[X1(:) Y1(:)]);       % evaluate eigenvectors on grid
idx = kmeans(V1, size(V1,2),'Replicates',10);       % kmeans clustering

%% plot partition
figure(3); clf; 
surf(X1,Y1,reshape(idx,ny1,nx1)); view(2); shading flat
axis equal; axis tight; xlabel('$x$'); ylabel('$y$'); colorbar

%% advect abd plot LCS
figure(4); clf; hold on; caxis([1 nc])
T = @(x) at_tf(flowmap(@bickleyjet,x,[t0 tf]));
for l = 2:nc
    I = find(idx==l); 
    S = [X1(I) Y1(I)]; 
    TS = T(S); TS(:,1) = mod(TS(:,1),20);
    scatter(TS(:,1),TS(:,2),10,'filled'); 
end
view(2); axis equal; axis([0 20 -3 3]); 
xlabel('lon [$^\circ$]'); ylabel('lat [$^\circ$]');


