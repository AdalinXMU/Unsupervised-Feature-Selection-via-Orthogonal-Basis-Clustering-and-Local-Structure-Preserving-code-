function V = update_V(W,X,E,mu,lambda1,U,V_old)
H = W'*X-E +lambda1/mu;
H = H* U;
% H= fillmissing(H,'previous');
inf_ind = isinf(H);
nan_ind = isnan(H);
temp = sum(inf_ind) + sum(nan_ind);
if temp> 0
    disp('nan apperance')
    V = V_old;
else
    [a,~,c] = svd(H,'econ');
    V = a * c';
end
