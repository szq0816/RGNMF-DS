function D_L21 = L21D(X,U,V)
%L21D 此处显示有关此函数的摘要
%   此处显示详细说明
D = [];
nsmp = size(X,2);
for i = 1:nsmp
    D(i) = 1/norm(X(:,i) - U * V(:,i));
end
    D_L21 = diag(D);
end

