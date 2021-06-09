function [W] = update_W_adaptive_inner_product(W,X,V,U,E,lambda1,mu,eta)

Q = V*U' + E - lambda1/mu;
p=1; % L_2p
[m,~] = size(W);
D = ones(m,m) - eye(m,m);
A =  (X * X') + 2*eta/mu*D;
W = pinv(A)*X*Q';


