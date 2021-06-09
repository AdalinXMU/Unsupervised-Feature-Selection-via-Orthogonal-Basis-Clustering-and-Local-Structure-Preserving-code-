function [idx,res] = Orthgonal_local_simliarity_without_simlarity_method(X,W,S_init,V,U,alpha,beta,eta,gamma)

NITER = 100;
A = (S_init+S_init')/2;
L_D = diag(sum(A));
L = L_D-A;
res = zeros(NITER,1);
res_one = norm(W'*X-V*U','fro')^2+eta*norm_21(W);
Y = W'*X;
size_n = size(S_init,1);
res_two = gamma*(trace(Y*L*Y' )+beta *norm(rand(size_n,size_n)-S_init,'fro')^2);
res_old = res_one + res_two ;
res(1) = res_old;
[m,~] = size(W);
D = eye(m);
p =1;
for iter=1:NITER
    Z = max(U,0);
    H = W'*X*U;
    [a,~,c] = svd(H,'econ');
    V = a * c';
    %     [W,D] = update_W_21norm(X,L,D,U,V,gamma,eta);
    temp_W = (X*X')+gamma*(X*L*X')+eta*D;
    W = pinv(temp_W)*X*U*V';
    tempD = 0.5*p * (sqrt(sum(W.^2,2)+eps)).^(p-2);
    D = diag(tempD);
    Y = W'*X;
    % S local S golbal and additive eps
%     S = new_update_S(S_init,Y,Y,beta);
%     A = (S+S')/2;
%     L_D = diag(sum(A));
%     L = L_D-A;
    H = X'*W*V + alpha*Z;
    [a,~,c] = svd(H,'econ');
    U = a * c';
    res_one = norm(W'*X-V*U','fro')^2+eta*norm_21(W);
    Y = W'*X;
    res_two = gamma*(trace(Y*L*Y' )+beta *norm(rand(size_n,size_n)-S_init,'fro')^2);
    res_new = res_one + res_two ;
    diff = res_old - res_new;
    res(iter+1) = res_new;
    if (iter > 1 && abs(diff)/(res_old) < 10^-4) ||(iter > 1 && abs(diff) < 10^-4);
        iter
        break
    else
        res_old = res_new;
    end
end

% plot(res)
[~,idx] = sort(sum(W.*W,2),'descend');

