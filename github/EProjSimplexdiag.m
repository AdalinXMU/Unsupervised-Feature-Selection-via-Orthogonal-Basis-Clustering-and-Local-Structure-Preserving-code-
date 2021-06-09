function [s, f] = EProjSimplexdiag(u, d)

%
%% Problem

%  min  1/2*s'*U*s-s'*d
%  s.t. s>=0, 1's=1


eta = min(u-d);
f = 1;
count=1;
while abs(f) > 10^-8
    v1 = 1./u*eta+d./u;
    posidx = v1>0;
    g = sum(1./u(posidx));
    f = sum(v1(posidx))-1;
    eta = eta - f/g;
    
    if count > 1000
        break;
    end;
    count=count+1;
end;
v1 = 1./u*eta+d./u;
s = max(v1,0);
