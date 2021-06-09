function [S,U] = ADA_ORTH_init(X, nClass)
%construct the affinity matrix
t = optSigma(X)^2;
S = constructW(X, struct('k',5, 'WeightMode', 'HeatKernel', 't', t));
diag_ele_arr = sum(S);
diag_ele_arr_t = diag_ele_arr.^(-1/2);
L = eye(size(X,1)) - diag(diag_ele_arr_t)* S *diag(diag_ele_arr_t);
L = (L + L')/2;
[eigvec, eigval] = eig(L);
[~, t1] = sort(diag(eigval), 'ascend');
eigvec = eigvec(:, t1(1:nClass));
eigvec = bsxfun(@rdivide, eigvec, sqrt(sum(eigvec.^2,2) + eps));

%init F and W
rand('twister',5489); %#ok
label = litekmeans(eigvec,nClass,'Replicates',10); % significantly!
U = rand(size(X,1),nClass);
for i = 1:size(X,1)
    U(i,label(i)) = 1;
end
U = U + 0.2;


end