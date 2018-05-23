addpath('../../src'); clear all

%% ABC flow 
t0 = 0; tf = 1; tspan = [t0,tf];
invCG = @(x) inv_CG_3d(@ABC,x,tspan);  % inverse Cauchy-Green tensor

%% nodes
n = 25; x = (0:n-1)/(n-1);
[p1,p2,p3] = meshgrid(x,x,x);
p = [p1(:),p2(:),p3(:)]; 
np = reshape(1:(n-1)^3,n-1,n-1,n-1); % periodic boundary conditions
np(n,:,:) = np(1,:,:); 
np(:,n,:) = np(:,1,:); 
np(:,:,n) = np(:,:,1);
pb = [(1:n^3)' np(:)];                   

%% triangulation
tri = delaunayTriangulation(p); 
t = tri.ConnectivityList;

%% assembly
deg = 1;  % quad degree
tic; G = inv_CG_quad_3d(p,t,invCG,deg); toc
tic; [D,M] = assemble_3d(p,t,pb,G); toc

%% eigenproblem
tic; [V,L] = eigs(D,M,10,'SM'); toc
[lam,ord] = sort(diag(L),'descend'); V = V(:,ord); 

%% plot spectrum
figure(1); clf; plot(lam,'*'); axis tight, axis square
xlabel('$k$'); ylabel('$\lambda_k$')

%% plot eigenvector
figure(2); clf; clear w; v = V(:,2);  theta = 0.2;
w(pb(:,1)) = v(pb(:,2)); 
plot_ev_3d(p1,p2,p3,w,theta);
axis equal; axis([0 1 0 1 0 1]); view(-145,25); caxis([-1,1])
camlight right; camlight left; lighting gouraud; colorbar
xlabel('$x$'); ylabel('$y$'); zlabel('$z$');

%% interpolate eigenvectors
nx1 = 70; x1 = linspace(0,1,nx1); 
[X1,Y1,Z1] = meshgrid(x1,x1,x1); 
W(pb(:,1),:) = V(pb(:,2),:);
V1 = eval_p1_3d(p1,p2,p3,W,[X1(:) Y1(:) Z1(:)]);      

%% partition by kmeans clustering
k = 5;
tic; idx = kmeans(V1(:,2:k),k,'Replicates',10); toc   

%% plot coherent sets
figure(3); clf;  hold on; %cols = colormap(jet(nc));
cols = ['m','b','y','g','r','c','k'];
h = hist(idx,k);[~,kmax] = max(h); 
for kk = 1:k,
    if kk~=kmax,
        plot_lcs(X1,Y1,Z1,idx,kk,cols(kk),0.7);
    end
end
axis equal; axis([0 1 0 1 0 1]); view(-145,25); 
camlight right; camlight left; lighting gouraud; 
xlabel('$x$'); ylabel('$y$'); zlabel('$z$');



