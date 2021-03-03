function Y = rowapply(f,x)

for i = 1:size(x,1),
    Y(:,:,i) = f(x(i,:));
end

