function [X,W,S,V,U,Z,E] = Init_convergence_adaptive_graph(X,class_num)

k = 5;
num = size(X,2);
dim = size(X,1);
d = ceil(dim/3);
d = min(d,num);

maxValue = max(max(X));
X = X/maxValue;
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

idx = kmeans(X',class_num);
U = zeros(num,class_num);

for i=1:num
    index = idx(i,1);
    U(i,index) = 1;
end
U = U + 0.2;
W = ones(dim,d);
V = W'*X*U;
E = W'*X - V*U';
Z = max(0,U);
