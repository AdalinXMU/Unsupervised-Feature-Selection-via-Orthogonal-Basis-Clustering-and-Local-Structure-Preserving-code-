clear;clc;
folder ='C:\Users\13263\Desktop\data';
addpath( genpath(folder) );
% dataset = 'ecoli_uni.mat';
dataset = 'Isolet1.mat';
disp(['dataset:',dataset]);
[X,Y] = extractXY(dataset);
X = double(X);
[nSmp,nDim] = size(X);
X = X';
class_num = length(unique(Y));
[X,W,S_init,V,U,Z,E] = Init_adaptive_graph(X,class_num);
gammaCandi= 10.^[-3:3];
betaCandi = 10.^[-3:3];
etaCandi =10.^[-3:3];
paramCell = fs_unsup_mymodel_build_param(gammaCandi,betaCandi,etaCandi);
disp('adaptive_graph_model ...');
t_start = clock;
feaSubsets = cell(length(paramCell), 1);
res_analysis = cell(length(paramCell), 1);
gamma = 0.01;
beta = 1;
eta = 0.1;
alpha = 10000;


[idx,res,S] = Orthgonal_local_simliarity_method(X,W,S_init,V,U,alpha,beta,eta,gamma);

res_d = [0;res];
res = [res;0];
plot(res-res_d)
