function S = new_update_S(S_init,U,Z,beta)
num = size(S_init,1);
dist = L2_distance_1(Z,U);
S = zeros(num);
islocal = 1;
index_P = S_init + S_init';
for i=1:num
    a0 = S_init(i,:);
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


