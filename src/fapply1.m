function Y = fapply1(f,X)

for i = 1:size(X,3),
    Y(:,:,i) = f(squeeze(X(:,:,i)));
end

