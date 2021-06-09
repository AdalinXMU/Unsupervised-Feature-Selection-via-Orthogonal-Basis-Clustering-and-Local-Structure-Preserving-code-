function [U,H] = update_U(W,X,E,Z,V,L,lambda1,lambda2,gamma,mu,U_old)

H0 = W'*X-E+lambda1/mu;
H1 = H0'*V;
H2 = Z  - gamma/mu*L*Z+ lambda2/mu;
H = (H1 + H2)*0.5;
inf_ind = isinf(H);
nan_ind = isnan(H);
temp = sum(inf_ind) + sum(nan_ind);
if temp> 0
    disp('inf apperance')
    U = U_old;
else
    [a,~,c] = svd(H,'econ');
    U = a * c';
end



