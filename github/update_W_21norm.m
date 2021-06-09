function [W,D] = update_W_21norm(X,L,D,U,V,gamma,eta)

INTER_W = 10;
% D = eye(size_n);
p=1; % L_2p
for i = 1:INTER_W
    temp_W = (X*X')+gamma*(X*L*X')+eta*D;
    W = pinv(temp_W)*X*U*V';
    tempD = 0.5*p * (sqrt(sum(W.^2,2)+eps)).^(p-2);
    D = diag(tempD);
    w1(i) = gamma*trace(W'*X*L*X'*W)+eta*sum(sqrt(sum(W.^2,2))); % log Tr(WXLXW)
    w2(i) = norm(W'*X-V*U','fro')^2;% gama*||W||_21
    WResult(i) = w1(i)+w2(i);
    if i > 1 && abs(WResult(i-1)-WResult(i)) < 0.000001
        break;
    end;
    
end;

end