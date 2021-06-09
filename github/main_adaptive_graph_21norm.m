clear;clc;
folder ='C:\Users\13263\Desktop\data';
addpath( genpath(folder) );
FeaNumCandi = [3,6,9,12,15,18]
% FeaNumCandi = [5,10,15,20,25,30,35,40,45,50]

% FeaNumCandi = [50,100,150,200,250,300]
prefix_mdcs = [];
nKmeans = 20;
% dataset = 'ORL_32x32';
% dataset = 'Yale_32x32';
% dataset = 'COIL20.mat';
% dataset = 'warpPIE10P';
% dataset = 'new_jaffe';
% dataset = 'new_UMIST';
dataset = 'lung_discrete';
% dataset = 'ecoli_uni.mat';
% dataset = 'dermatology_uni';
% dataset = 'ionosphere_uni';
% dataset = 'colon';
% dataset = 'binalpha_uni';
% dataset = 'MSRA25_uni';
% dataset = 'Isolet1.mat';
disp(['dataset:',dataset]);
[X,Y] = extractXY(dataset);
X = double(X);
[nSmp,nDim] = size(X);
X = X';
class_num = length(unique(Y));

[X,W,S_init,V,U,Z,E]= Init_adaptive_graph3(X,class_num);
% betaCandi = 10;
% etaCandi =100;
% gammaCandi =10;
% gammaCandi= 10.^[-5:5];
gammaCandi= 10.^[-3:3];
betaCandi = 10.^[-3:3];
etaCandi =10.^[-3:3];
paramCell = fs_unsup_mymodel_build_param(gammaCandi,betaCandi,etaCandi);
disp('adaptive_graph_model ...');
t_start = clock;
feaSubsets = cell(length(paramCell), 1);
res_analysis = cell(length(paramCell), 1);
S_analysis = cell(length(paramCell), 1);
% alpha = 100;
for i1 = 1:length(paramCell)
    fprintf(['adaptive_graph model parameter search %d out of %d...\n'], i1, length(paramCell));
    param = paramCell{i1};
    gamma = param.gamma;
    beta = param.beta;
    eta = param.eta;
    alpha = max([gamma,beta,eta])*10
    [idx,res,S_iter] = Orthgonal_local_simliarity_method(X,W,S_init,V,U,alpha,beta,eta,gamma);
    feaSubsets{i1,1} = idx;
    res_analysis{i1,1} = res;
    S_analysis{i1,1} = S_iter;
end
t_end = clock;
t1 = etime(t_end,t_start);
disp(['exe time: ',num2str(t1)]);
disp('evaluation....');
t_start = clock;
res_aio = cell(length(paramCell), length(FeaNumCandi));
for i2 = 1:length(FeaNumCandi)
    parfor i1 = 1:length(paramCell)
        fprintf('mymodel parameter evaluation %d outof %d  ... %d out of %d...\n', i2, length(FeaNumCandi), i1, length(paramCell));
        idx = feaSubsets{i1,1};
        res_aio{i1, i2} = evalUnSupFS(X', Y, idx(1:FeaNumCandi(i2)), struct('nKm', nKmeans));
    end
end
[res_gs, res_gs_ps] = grid_search_fs(res_aio);
res_gs.feaset = FeaNumCandi;
t_end = clock;
t2 = etime(t_end,t_start);
disp(['exe time: ',num2str(t2)]);
res_gs.time = t1;
res_gs.time2 = t2;

save(fullfile(prefix_mdcs, [dataset, '_best_result_ada_orth_21norm.mat']), 'FeaNumCandi','res_gs','res_aio', 'res_gs_ps');
[max_acc,ind] = max(res_gs.mean_acc);
std_acc = res_gs.std_acc(ind);
[max_nmi,index] = max(res_gs.mean_nmi_sqrt);
std_nmi = res_gs.std_nmi_sqrt(index);
result = [max_acc,std_acc,max_nmi,std_nmi].*100;
result = roundn(result,-2)
% sound(sin(2*pi*25*(1:4000)/100));
sound(sin(2*pi*25*(1:4000)/100));
% res_temp =res_analysis{1,1};
% S_temp = S_analysis{1,1};
% HeatMap(S_temp)
% plot(res_temp(1:4))
% plot(res_temp(3:20))
% plot(res_analysis{1,1})