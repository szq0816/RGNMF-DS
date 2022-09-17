function similarity = ConstructS(data, labs, n_labeled)
%CONSTRUCTD 此处显示有关此函数的摘要
%   此处显示详细说明

% n_unlabeled = length(labs) - n_labeled;

labeled = labs(1:n_labeled);

% D = labeled' * (1 ./labeled);
S = eye(length(labs));

for i = 1:n_labeled
    for j = 1:n_labeled
        if labeled(i)==labeled(j)
           S(i,j) = 1;
        else
           S(i,j)= 0;
        end
    end

end

options.NeighborMode = 'KNN';
options.k = 5;
options.WeightMode = 'HeatKernel';
options.t = 100;

SS = constructW(data,options);
SS(1:n_labeled, 1:n_labeled) = 0;

similarity = S+SS;

end

