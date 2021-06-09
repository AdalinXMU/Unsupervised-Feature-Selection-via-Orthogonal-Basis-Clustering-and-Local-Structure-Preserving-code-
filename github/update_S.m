function S = update_S(A,U,Z,beta)
num = size(A,1);
dist = L2_distance_1(Z',U');
S = zeros(num);
islocal = 1;
index_P = A + A';
for i=1:num
    a0 = A(i,:);
    a1 = index_P(i,:);
    if islocal == 1
        idxa0 = find(a1>0);
    else
        idxa0 = 1:num;
    end;
    ai = a0(idxa0);
    di = dist(i,idxa0);
    ad = ai-1/(2*beta)*di;
    S(i,idxa0) = EProjSimplex_new(ad,1);
end;


