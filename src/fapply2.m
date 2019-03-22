function Y = fapply2(f,X,W)

for j = 1:size(X,3)
    Y(:,:,j) = f(squeeze(X(:,:,j)),squeeze(W(:,:,j)));
end

