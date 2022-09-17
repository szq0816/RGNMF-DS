function [U,V,obj] = rssNMF(X,para,Q,D, initU,initV, initS)
% X in d \times N

alpha=para.alpha;  % the hyper papramter lambda
k=para.k;
maxiter=para.maxiter;
beta=para.beta;

d=size(X,1);
n=size(X,2);
L=D-Q;

% init
% U=rand(d,k);
% V=rand(k,n);

U = initU;
V = initV;
S = initS;


obj(1)=sum(sum((X-U*V-S).^2))+alpha*norm(S,1)+beta*trace(V*L*V');
for iter=1:maxiter
    
    U=U.*((X * V')./(U*V *V' + S*V' + eps));
    V=V.*((U'*X + beta * V * Q)./(U'*U*V + U' * S + beta * V * D + eps));
    S = soft(X-U*V, 0.5 * alpha);

    % normization
    V=V.*(repmat(sum(U,1)',1,n));
    U=U./(repmat(sum(U,1),d,1));
%     S=S./(repmat(sum(S,1),d,1));
    
    obj(iter+1)=sum(sum((X-U*V-S).^2))+alpha*norm(S,1)+beta*trace(V*L*V');
    disp(['the ', num2str(iter), ' obj is ', num2str(obj(iter))]);
    if ((abs(obj(iter)-obj(iter+1)) / obj(iter))<10^-3)
        break;
    end
end
end

