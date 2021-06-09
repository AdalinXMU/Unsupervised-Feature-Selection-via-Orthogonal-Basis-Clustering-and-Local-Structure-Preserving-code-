function [V,E] = Init_V_E(X,W,U)

H = U'*X'*W;

[a,~,c] = svd(H,'econ');
% V = a * c';
V = c*a';
E = W'*X - V*U';