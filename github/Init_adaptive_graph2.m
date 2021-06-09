function [X,W,S,V,U,Z,E] = Init_adaptive_graph2(X,nClass)
k = 5;
gamma = 10;
num = size(X,2);
dim = size(X,1);
% d = ceil(dim/2);
% d = min(d,num);
d = nClass;
maxValue = max(max(X));
X = X/maxValue;
[~,U] = ADA_ORTH_init(X', nClass);
X0 = X';
mX0 = mean(X0);
X1 = X0 - ones(num,1)*mX0;
scal = 1./sqrt(sum(X1.*X1)+eps);
scalMat = sparse(diag(scal));
X = X1*scalMat;
X = X';

distX = L2_distance_1(X,X);
[distX1, idx] = sort(distX,2);
S = zeros(num);
rr = zeros(num,1);
for i = 1:num
    di = distX1(i,2:k+2);
    rr(i) = 0.5*(k*di(k+1)-sum(di(1:k)));
    id = idx(i,2:k+2);
    S(i,id) = (di(k+1)-di)/(k*di(k+1)-sum(di(1:k))+eps);
end;
r = mean(rr);
lambda = r;

A0 = (S+S')/2;
D0 = diag(sum(A0));
L0 = D0 - A0;
[F, ~, evs]=eig1(L0, nClass, 0);

[W] = InterationW(L0,X,gamma,dim,d);
NITER = 50;
for iter = 1:NITER
    distf = L2_distance_1(F',F');
    distx = L2_distance_1(W'*X,W'*X);
    if iter>5
        [~, idx] = sort(distx,2);
    end;
    S = zeros(num);
    for i=1:num
        idxa0 = idx(i,2:k+1);
        dfi = distf(i,idxa0);
        dxi = distx(i,idxa0);
        ad = -(dxi+lambda*dfi)/(2*r);
        S(i,idxa0) = EProjSimplex_new(ad);
    end;
    
    S = (S+S')/2;
    D = diag(sum(S));
    L = D-S;
    
    [W] = InterationW(L,X,gamma,dim,d);
    F_old = F;
    [F, ~, ev]=eig1(L, nClass, 0);
    evs(:,iter+1) = ev;
    
    fn1 = sum(ev(1:nClass));
    fn2 = sum(ev(1:nClass+1));
    if fn1 > 0.000000001
        lambda = 2*lambda;
    elseif fn2 < 0.00000000001
        lambda = lambda/2;  F = F_old;
    elseif iter>1
        break;
    end;
end;

% W = eye(dim,d); 
V = W'*X*U;
E = W'*X - V*U';
Z = max(0,U);
