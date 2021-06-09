function [W] = update_W_adaptive_similarity(W,X,V,U,E,lambda1,mu,eta)
% W = max(W,10^7);
% Q = V*U' + E + lambda1/mu;
Q = V*U' + E - lambda1/mu;
INTER_W = 1;
p=1; % L_2p

for i = 1:INTER_W
    tempD = 0.5*p * (sqrt(sum(W.^2,2)+eps)).^(p-2);
    D = diag(tempD);
    A =  (X * X') + 2*eta/mu*D;
    W = pinv(A)*X*Q';
    w1(i) = mu/2*norm(W'*X - Q,'fro')^2; % log Tr(WXLXW)
    w2(i) = eta*sum(sqrt(sum(W.^2,2)));% gama*||W||_21
    WResult(i) = w1(i)+w2(i);
    if i > 1 && abs(WResult(i-1)-WResult(i)) < 0.000001
        break;
    end;
end;

end