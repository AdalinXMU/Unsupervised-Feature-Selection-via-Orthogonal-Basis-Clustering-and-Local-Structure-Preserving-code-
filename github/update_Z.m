function Z = update_Z(L,U,gamma,mu,lambda2)

P = U - lambda2/mu - gamma*L*U/mu;

Z = max(P,0);

