function Dissimilarity = ConstructD(labs,n_labeled)
%CONSTRUCTD 此处显示有关此函数的摘要
%   此处显示详细说明

% n_unlabeled = length(labs) - n_labeled;

labeled = labs(1:n_labeled);

% D = labeled' * (1 ./labeled);
Dissimilarity = zeros(length(labs));

for i = 1:n_labeled
    for j = 1:n_labeled
        if labeled(i)==labeled(j)
           Dissimilarity(i,j) = 0;
        else
           Dissimilarity(i,j)= 1;
        end
    end

end


end

