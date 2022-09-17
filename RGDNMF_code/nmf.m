function [W, H] = nmf(X, K, MAXITER,initU,initV)
%Euclidean distance
[m,n]= size(X);

if ~isempty(initU)
    W = initU;
    H = initV;
    
else

    rand('seed',0)
    W = abs(rand(m, K));
    H = abs(rand(K, n));
end

dnorm0 = norm(X - W * H, 'fro');
for i=1:MAXITER
    W = W .* (X*H')./(W*H*H'+eps);
    H = H .* (W' * X)./(W'*W*H+eps) ; 
    dnorm = norm(X - W * H, 'fro');
    err = abs((dnorm0-dnorm) / max(1,dnorm0));
    if err <= 0.05
        fprintf('out: %d',i);
        break;
    end
    dnorm0 = dnorm;
end