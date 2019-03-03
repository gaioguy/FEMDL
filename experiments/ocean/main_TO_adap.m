addpath('../../src','../../src/2d'); clear all; clc; colormap jet

%% ocean vector field and flow map
% The altimeter data are produced by SSALTO/DUACS and distributed by 
% AVISO, (http://www.aviso.oceanobs.com/duacs).
load('Ocean_geostrophic_velocity.mat','lon','lat','UT','VT','time');
U = griddedInterpolant({lon,lat,time},permute(UT,[2,1,3]),'cubic','none');
V = griddedInterpolant({lon,lat,time},permute(VT,[2,1,3]),'cubic','none');

%% flow map
vf = @(t,x) ocean(t,x,U,V);
t0 = time(1); tf = t0+90; nt = 2; tspan = linspace(t0,tf,nt);
T = @(x) flow_map(vf,x,tspan);

%%  nodes
xmin = -4; xmax = 6; ymin = -34; ymax = -28; 
nx = 150; ny = 0.6*nx; n = nx*ny; 
x1 = linspace(xmin,xmax,nx); y1 = linspace(ymin,ymax,ny);
[X,Y] = meshgrid(x1,y1); 
p0 = [X(:) Y(:)];          % nodes
pb = [1:n; 1:n]';                                   % no periodic boundary

%% triangulation
tri = alphaShape(p0(:,1),p0(:,2));
b = unique(boundaryFacets(tri));                    % b = boundary nodes

%% time integration
tic; P = T(p0); toc
for k = 1:nt, p{k} = [P(:,k) P(:,k+nt)]; end;

%% assembly (for missing data case)
pm = 1; tic                                    % percentage of nodes to remove
D = sparse(n,n); M = sparse(n,n); 
for k = 1:nt
    r = randperm(n,floor(pm*n))'; 
    tri = alphaShape(p{k}(r,1),p{k}(r,2),0.4);
    tr = alphaTriangulation(tri); 
    t = [r(tr(:,1)), r(tr(:,2)), r(tr(:,3))];
    I = kron([1 0 1],ones(size(t,1),1));         % 2 x 2 identity matrix
    [Dt,Mt{k}] = assemble(p{k},t,pb,I);
    D = D + Dt/nt; M = M + Mt{k}/nt; 
end
toc

%% Dirichlet boundary condition
D(b,:) = 0; D(:,b) = 0; M(b,:) = 0; M(:,b) = 0; 
D(b,b) = speye(length(b),length(b));

%% remove all zero rows and columns
S = sum(abs(D)); I = find(abs(S)>eps); 
D = D(I,I); M = M(I,I); pI = p{1}(I,:); pbI = [1:size(pI,1); 1:size(pI,1)]';

%% solve eigenproblem
tic; [V,L] = eigs(D,M,10,'SM'); toc
[lam,ord] = sort(diag(L),'descend'); V = V(:,ord);

%% plot spectrum
figure(1); clf; plot(lam,'*'); axis tight, axis square
xlabel('$k$'); ylabel('$\lambda_k$')

%% plot eigenvectors
figure(2), clf, 
tri = alphaShape(pI(:,1),pI(:,2)); 
tI = alphaTriangulation(tri);
plotf(pI,tI,pbI,V(:,1),0); 
xlabel('lon [$^\circ$]'); ylabel('lat [$^\circ$]'); 

%% compute partition
nx1 = 400; ny1 = 0.6*nx1; 
x1 = linspace(xmin,xmax,nx1); y1 = linspace(ymin,ymax,ny1);
[X1,Y1] = meshgrid(x1,y1); 
nc = 3;
tic; V1 = eval_p1(pI,V(:,1:nc),[X1(:) Y1(:)]); toc        % evaluate eigenvectors on grid
tic; idx = kmeans(V1, nc+1, 'Replicates', 10); toc          % kmeans clustering

%% plot partition
figure(3); clf; surf(X1,Y1,reshape(idx,ny1,nx1)); view(2); shading flat
axis equal; axis tight; xlabel('lon [$^\circ$]'); ylabel('lat [$^\circ$]');  
load cmap7; colormap(cmap); 

%% advect abd plot LCS
figure(4); clf; hold on; colormap(cmap); caxis([1 nc+1])
T = @(x) flow_map(vf,x,tspan);
for l = 2:nc+1
    I = find(idx==l); 
    S = [X1(I) Y1(I)];
    TS = T(S); TS = TS(:,[2,4]);
    scatter3(TS(:,1),TS(:,2),zeros(length(I),1),10,cmap(l,:),'filled'); 
end
view(2); axis equal; axis([-8 2 -33 -27]); colorbar
xlabel('lon [$^\circ$]'); ylabel('lat [$^\circ$]');
