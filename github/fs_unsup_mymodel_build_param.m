function paramCell = fs_unsup_mymodel_build_param(gammaCandi,betaCandi,etaCandi)
n1 = length(gammaCandi);
n2 = length(betaCandi);
n3 = length(etaCandi);

nP = n1*n2*n3;
paramCell = cell(nP, 1);
idx = 0;
for i1 = 1:n1
    for i2 = 1:n2
        for i3 = 1:n3
            param = [];
            param.gamma = gammaCandi(i1);
            param.beta = betaCandi(i2);
            param.eta = etaCandi(i3);
            idx = idx + 1;
            paramCell{idx} = param;
        end
    end
end