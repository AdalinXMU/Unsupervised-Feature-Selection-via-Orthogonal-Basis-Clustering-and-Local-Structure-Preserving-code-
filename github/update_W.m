function [W,P] = update_W(X,V,U,E,S,lambda1,mu,beta,eta)

[size_n,size_m] = size(X);
% P = zeros(size_n,size_n);
S_square = S.*S;
S_square = (S_square+ S_square')/2;
D_square = diag(sum(S_square));
L_square = D_square - S_square;
P = X*L_square*X';

Q = V*U' + E + lambda1/mu;
A = beta * P + mu/2 * (X * X') + eta*(ones(size_n) - eye(size_n));
W = mu/2 * pinv(A)*X*Q';
end