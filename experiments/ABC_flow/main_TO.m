addpath('../../src'); clear all

%% ABC inverse flow map
t0 = 0; tf = 1; tspan = [tf,t0]; 
T1 = @(x) flow_map_3d(@ABC,x,tspan);

%% nodes
n = 25; x = (0:n-1)/(n-1);
[p1,p2,p3] = meshgrid(x,x,x);
p = [p1(:),p2(:),p3(:)]; 
np = reshape(1:(n-1)^3,n-1,n-1,n-1); % periodic boundary conditions
np(n,:,:) = np(1,:,:); 
np(:,n,:) = np(:,1,:); 
np(:,:,n) = np(:,:,1);
pb = [(1:n*n*n)' np(:)];        % nodes pb(:,1) are identified with pb(:,2)

%% triangulation
tri = delaunayTriangulation(p); 
t = tri.ConnectivityList;

%% integration and assembly
tic; P = T1(p); T1p = P(:,[2,4,6]); toc
A = kron([1 0 0 1 0 1],ones(size(t,1),1));  % 3 x 3 idenitiy matrix
tic; [D,M] = assemble_3d(p,t,pb,A); toc

%% compute transfer operator and dynamic Laplacian
Phi = get_Phi_3d(p,t,pb);
tic; Alpha = get_alpha_3d(tri,pb,T1p,Phi); toc
DL = 0.5*(D + Alpha'*D*Alpha);

%% eigenproblem
tic; [V,L] = eigs(DL,M,10,'SM'); toc
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
tic; V1 = eval_p1_3d(p1,p2,p3,W,[X1(:) Y1(:) Z1(:)]); toc      

%% partition by kmeans clustering
nc = 5;
tic; idx = kmeans(V1(:,1:nc-1),nc,'Replicates',10); toc   

%% plot coherent sets
figure(3); clf;  hold on; %cols = colormap(jet(nc));
cols = ['y','r','b','g','m','c','k'];
h = hist(idx,nc);[~,kmax] = max(h); 
for k = 1:nc,
    if k~=kmax,
        plot_lcs(X1,Y1,Z1,idx,k,cols(k),0.7);
    end
end
axis equal; axis([0 1 0 1 0 1]); view(-145,25); 
camlight right; camlight left; lighting gouraud; 
xlabel('$x$'); ylabel('$y$'); zlabel('$z$');
