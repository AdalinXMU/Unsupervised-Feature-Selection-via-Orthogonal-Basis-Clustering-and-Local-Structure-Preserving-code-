function  [L0, alpha] = determined_alpha_initilize_S(X,num,k)

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
alpha = mean(rr);

A0 = (S+S')/2;
D0 = diag(sum(A0));
L0 = D0 - A0;  