function [H,temp] = update_test(W,X,E,Z,V,L,lambda1,lambda2,gamma,mu)

H0 = W'*X-E-lambda1/mu;
H1 = H0'*V;

H2 = Z  - gamma/mu*L*Z;
temp = L*Z;
H2 = H2+ lambda2/mu;
H = H1 + H2;







