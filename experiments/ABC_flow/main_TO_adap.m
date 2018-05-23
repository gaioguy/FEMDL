addpath('../../src/3d'); clear all

%% ABC flow map
tspan = [0,1]; nt = 2;
T = @(x) flow_map_3d(@ABC,x,tspan);

%% data points
n = 25; x = (0:n-2)/(n-1);
[px,py,pz] = meshgrid(x,x,x);
P0 = [px(:),py(:),pz(:)]; 

%% time integration
tic; P = T(P0); toc
P1 = [P(:,2) P(:,4) P(:,6)]; 

%% triangulations
[p0,t0,pb0] = delaunay_T3(P0, [0 1 0 1 0 1]);
[p1,t1,pb1] = delaunay_T3(P1, [0 1 0 1 0 1]);

%% assembly
A = kron([1 0 0 1 0 1],ones(size(t0,1),1));     % 3 x 3 identity matrix
tic; [D0,M] = assemble_3d(p0,t0,pb0,A); toc
A = kron([1 0 0 1 0 1],ones(size(t1,1),1));     % 3 x 3 identity matrix
tic; [D1,M1] = assemble_3d(p1,t1,pb1,A); toc
D = (D0 + D1)/2; 

%% eigenproblem
tic; [V,L] = eigs(D,M,10,'SM'); toc
[lam,ord] = sort(diag(L),'descend'); V = V(:,ord);

%% plot spectrum
figure(1); clf; plot(lam,'*'); axis tight, axis square
xlabel('$k$'); ylabel('$\lambda_k$')

%% plot eigenvector
figure(2); clf; clear w; v = -V(:,4);  theta = 0.4;
w(pb0(:,1)) = v(pb0(:,2)); 
plot_ev_3d(px,py,pz,v,theta);
axis equal; axis([0 1 0 1 0 1]); view(-145,25); caxis([-1,1])
camlight right; camlight left; lighting gouraud; colorbar
xlabel('$x$'); ylabel('$y$'); zlabel('$z$');

%% interpolate eigenvectors
nx1 = 70; x1 = linspace(0,1,nx1); 
[X1,Y1,Z1] = meshgrid(x1,x1,x1); 
V1 = eval_p1_3d(px,py,pz,V,[X1(:) Y1(:) Z1(:)]);      

%% partition by kmeans clustering
nc = 5;
tic; idx = kmeans(V1(:,1:nc),nc,'Replicates',2); toc   

%% plot coherent sets
figure(3); clf;  hold on; % cols = colormap(jet(nc));
cols = ['b','r','g','y','m','c','c'];
h = hist(idx,nc);[~,kmax] = max(h); 
for k = 1:nc,
    if k~=kmax,
        plot_lcs(X1,Y1,Z1,idx,k,cols(k),0.7);
    end
end
axis equal; axis([0 1 0 1 0 1]); view(-145,25); 
camlight right; camlight left; lighting gouraud; 
xlabel('$x$'); ylabel('$y$'); zlabel('$z$');

26.122