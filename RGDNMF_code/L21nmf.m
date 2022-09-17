function [U,V,obj]=L21nmf(X,para, initU,initV)
% X in d \times N

lambda=para.lambda;  % the hyper papramter lambda
k=para.k;
maxiter=para.maxiter;
mu=para.mu;
alpha = para.alpha;

d=size(X,1);
n=size(X,2);

% init
% U=rand(d,k);
% V=rand(k,n);

U = initU;
V = initV;

D_L21 = L21D(X,U,V);

obj(1)=sum(sqrt(sum((X-U*V).^2)));
for iter=1:maxiter
    D_L21 = L21D(X,U,V);
    
    U=U.*((X * D_L21 * V')./(U*V * D_L21 *V'+eps));
    V=V.*((U'*X * D_L21)./(U'*U*V * D_L21 + eps));
    Z=X-U*V;

    % normization
    V=V.*(repmat(sum(U,1)',1,n));
    U=U./(repmat(sum(U,1),d,1));
    obj(iter+1)=sum(sqrt(sum(Z).^2));
    disp(['the ', num2str(iter), ' obj is ', num2str(obj(iter))]);
    if ((abs(obj(iter)-obj(iter+1)) / obj(iter))<10^-3)
        break;
    end
end