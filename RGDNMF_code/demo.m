clc
clear;
close all;
% 'Test_1_Zeisel_big',
dataset = {'Test_Kolod','Test_Darmanis','Test_Ramskold', 'Test_human', 'Test_islet', ...
    'Test_Zheng','GSE57249','GSE64016','GSE70657','Armstrong','bhattacharjee','lapointe',...
    'pomeroy','GSE81252','GSE75688','GSE85908','Quake_10x_Limb_Muscle','Test_1_Zeisel_big'};

initdataset = {'init_Test_Kolod','init_Test_Darmanis','init_Test_Ramskold', 'init_Test_human', 'init_Test_islet','init_Test_Zheng', ...
    'init_GSE57249','init_GSE64016','init_GSE70657','init_Armstrong','init_bhattacharjee','init_lapointe','init_pomeroy',...
    'init_GSE81252','init_GSE75688','init_GSE85908','init_Quake_10x_Limb_Muscle','init_Test_1_Zeisel_big'};
addpath('../data')
addpath('./initdata')
ndata = length(dataset);
collect_result_NMI = [];
collect_result_ARI = [];
collect_result_ACC = [];
collect_result_Purity = [];
data_num = [1,3,4,5,7,8,9, 14,15,16,17,18];
% data_num = [1:ndata];
result = [];
rng(5)
for i = 1:12
    
    load(dataset{data_num(i)});
    load(initdataset{data_num(i)});
%     initU = initU;
%     initV = initV;
%     rand_orad = randperm(length(true_labs));
%     rand_orad = rand_orad;
    X = in_X(rand_orad,:);
%     Xt = normalize(X');
    Xt = NormalizeFea(X);
    X = Xt;
    labs = true_labs(rand_orad);
    n_labeed = floor(0.2 * length(labs));
    nfeature = floor(sqrt(size(X,1)));
% X = in_X;
% labs = true_labs;

    nnClass = length(unique(true_labs)); 

    NMI = [];
    ARI = [];
    ACC = [];
    Purity= [];
    for j = 1:1
        fprintf('time: %d',j);
 
       
    %% nmf
    
%     [U, V] = nmf(X', floor(size(X,1)/2), 500,initU,initV);
%     [project_labs, center] = litekmeans(V', nnClass, 'Start', idx(1:nnClass));
% %     
% %     [U0, V0] = nmf(X', floor(size(X,1)/2), 500,initU,initV);
   %% L21nmf
    
%     [U, V, obj] = L21nmf(X', para, initU,initV);
%     [project_labs, center] = litekmeans(V', nnClass, 'Start', idx(1:nnClass));
%     
% %     [U0, V0] = nmf(X', floor(size(X,1)/2), 500,initU,initV);

    %% GNMF
%     addpath('../Tools')
%     options = [];
%     options.NeighborMode = 'KNN';
%     options.k = 5;
%     options.WeightMode = 'HeatKernel';
%     options.t = 10;
% 
%     W = constructW(X,options);%nnSWM1(X,X,5,1,0);
% 
%     options = [];
%     options.error = 0.005;
%     options.maxIter = 500;
%     options.nRepeat = 1;
%     options.minIter = 20;
%     options.meanFitRatio = 0.5;
%     options.alpha = 10;
% %     U = [];
% %     V = [];
%     
% %     k = floor(sqrt(nsmp));
% %     [U, V] = nmf(X,nfeature,200);
%     U = initU;
%     V = initV;
% %     
% %     [Ze, He] = GNMF_Multi(X', nfeature, W, options,U, V');
% %     [project_labs, center] = litekmeans(He, nnClass, 'Start', idx(1:nnClass));
%    
%     [U2, V2] = GNMF_Multi(X', nfeature, W, options,U, V');
    

    %% rssnmf
% %     para = [];
% %     para.alpha = 8;
% %     para.k = nfeature;
% %     para.maxiter = 200;
% %     para.beta = 150;
% 
% 
%     % Dissimilarity and Similarity matrix
% % 
% %     D = ConstructD(labs, n_labeed); %dissimilarity matrix
% %     S = ConstructS(X, labs, n_labeed); %similarity matrix
% %     A = diag(sum(S, 2)); %diagonal matrix
% %    
%     % Graph matrix
%     options = [];
%     options.NeighborMode = 'KNN';
%     options.k = 5;
%     options.WeightMode = 'HeatKernel';
%     options.t = 1;
% 
%     W = constructW(X,options);%nnSWM1(X,X,5,1,0);
%     D = diag(sum(W, 2));
% 
%     
%     % init
%     [d,n] = size(X');
% %     initU = rand(d,nfeature);
% %     initV = rand(nfeature,n);
%     initS = rand(d,n);
%     
% 
%     [U,V,obj] = rssNMF(X', para, W, D, initU,initV,initS);
% %     
%       [project_labs, center] = litekmeans(V', nnClass, 'Start', idx(1:nnClass));
%     
% 
% 


     %% RGNMF-DS
%     para = [];
%     para.lambda = 0;
%     para.k = nfeature;
%     para.maxiter = 100;
%     para.mu = 0;
%     para.alpha = 0;

    % Dissimilarity and Similarity matrix

    D = ConstructD(labs, n_labeed); %dissimilarity matrix
    S = ConstructS(X, labs, n_labeed); %similarity matrix
    A = diag(sum(S, 2)); %diagonal matrix
   
    % Graph matrix
    options = [];
    options.NeighborMode = 'KNN';
    options.k = 5;
    options.WeightMode = 'HeatKernel';
    options.t = 10;

    GW = constructW(X,options);%nnSWM1(X,X,5,1,0);
    GD = diag(sum(GW, 2));
    
    [U,V,obj] = RGNMF_DS(X', D, para, A, S, GW, GD, initU, initV);
    
      [project_labs, center] = litekmeans(V', nnClass, 'Start', idx(1:nnClass));

    %% 评估
    NMI(j) = Cal_NMI(labs, project_labs);
    ARI(j) = Cal_ARI(labs, project_labs);
    ACC(j) = ACC_ClusteringMeasure(labs, project_labs);
    [~,Purity0,~,~] = Fmeasure(labs', project_labs');
    Purity(j) = Purity0;
    NMIj = NMI(j);
    ARIj = ARI(j);
    ACCj = ACC(j);
    Purityj = Purity(j);
%     eval(['save', ' ', './initdata/init_', dataset{data_num(i)},'_' num2str(j), '.mat', ' ', 'initU initV rand_orad idx NMIj ARIj para'])
%     save ./initdata/init.mat initU initV rand_orad;
%     eval(['save', ' ', './obj_', dataset{data_num(i)},'_' num2str(j), '.mat', ' ', 'obj'])
%     save weight.mat X V D S GW labs;
%     eval(['save', ' ', '../plot/C_NMF_Graph_', dataset{data_num(i)},'_' num2str(j), '.mat', ' ', 'obj3'])
    end
    collect_result_NMI(i) = mean(NMI,2);
    collect_result_ARI(i) = mean(ARI,2);
    collect_result_ACC(i) = mean(ACC,2);
    collect_result_Purity(i) = mean(Purity,2);    

end
    result =[collect_result_ACC;collect_result_ARI;collect_result_NMI;collect_result_Purity]





