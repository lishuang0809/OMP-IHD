% This code is used to generate Table 2.
% Coded by Shuang Li, Jan. 2025.

clc; clear;

J = 100;            % # of samples
kk = [10 30 35 58 60 100];   % sparsity
NN = [64 64 128 128 256 256];
LL = 2*NN;     % The dictionary is of size N x L

%alg_all = 1:10;
rho_snr = 60;   %in dB, used to define a successful recovery.

Rate_OMP = zeros(1,length(kk));
time_OMP = zeros(1,length(kk));
Rate_CoSaMP = zeros(1,length(kk));
time_CoSaMP = zeros(1,length(kk));
Rate_BP = zeros(1,length(kk));
time_BP = zeros(1,length(kk));
Rate_Yin_IRLS = zeros(1,length(kk));
time_Yin_IRLS = zeros(1,length(kk));
Rate_FISTA = zeros(1,length(kk));
time_FISTA = zeros(1,length(kk));
Rate_GL2 = zeros(1,length(kk));
time_GL2 = zeros(1,length(kk));
Rate_GL1 = zeros(1,length(kk));
time_GL1 = zeros(1,length(kk));
Rate_FGL1 = zeros(1,length(kk));
time_FGL1 = zeros(1,length(kk));
Rate_GLQ = zeros(1,length(kk));
time_GLQ = zeros(1,length(kk));
Rate_FGLQ = zeros(1,length(kk));
time_FGLQ = zeros(1,length(kk));


%condi=10000;  [Q,S_opt,Q_condi]=gen_data_condQ(N,L,J,kappa,condi);

for i = 1:length(kk)
    k = kk(i);
    N = NN(i);
    L = LL(i);
    [X, Q, S_opt] = data_generation(N, L, k, J);  % X = Q*S_opt, ||X(:,j)||_2=1

    % ======= OMP_C ========
    [S_OMP,time_OMP(i)] = finalized_Master_4_Lis(X,Q,k,1);
    [Rate_OMP(i),~,~] = ComputingRecoveryRate(S_OMP,S_opt,Q,rho_snr);

    % ======= CoSaMP_C ========
    [S_CoSaMP,time_CoSaMP(i)] = finalized_Master_4_Lis(X,Q,k,2);
    [Rate_CoSaMP(i),~,~] = ComputingRecoveryRate(S_CoSaMP,S_opt,Q,rho_snr);

    % ======= BP_C ========
    [S_BP,time_BP(i)] = finalized_Master_4_Lis(X,Q,k,3);
    [Rate_BP(i),~,~] = ComputingRecoveryRate(S_BP,S_opt,Q,rho_snr);

    % ======= Yin_IRLS ========
    [S_Yin_IRLS,time_Yin_IRLS(i)] = finalized_Master_4_Lis(X,Q,k,4);
    [Rate_Yin_IRLS(i),~,~] = ComputingRecoveryRate(S_Yin_IRLS,S_opt,Q,rho_snr);

    % ======= FISTA ========
    [S_FISTA,time_FISTA(i)] = finalized_Master_4_Lis(X,Q,k,5);
    [Rate_FISTA(i),~,~] = ComputingRecoveryRate(S_FISTA,S_opt,Q,rho_snr);


    % ======= Alg_{GL2} ========
    [S_GL2,time_GL2(i)] = finalized_Master_4_Lis(X,Q,k,6);
    [Rate_GL2(i),~,~] = ComputingRecoveryRate(S_GL2,S_opt,Q,rho_snr);

    %======= Alg_{GL1} =======
    [S_GL1,time_GL1(i)] = finalized_Master_4_Lis(X,Q,k,7);
    [Rate_GL1(i),~,~] = ComputingRecoveryRate(S_GL1,S_opt,Q,rho_snr);

    %======= Alg_{GL1}^F =======
    [S_FGL1,time_FGL1(i)] = finalized_Master_4_Lis(X,Q,k,8);
    [Rate_FGL1(i),~,~] = ComputingRecoveryRate(S_FGL1,S_opt,Q,rho_snr);

    % ======= Alg_{GLQ} ========
    [S_GLQ,time_GLQ(i)] = finalized_Master_4_Lis(X,Q,k,9);
    [Rate_GLQ(i),~,~] = ComputingRecoveryRate(S_GLQ,S_opt,Q,rho_snr);

    % ======= Alg_{GLQ}^F ========
    [S_FGLQ,time_FGLQ(i)] = finalized_Master_4_Lis(X,Q,k,10);
    [Rate_FGLQ(i),~,~] = ComputingRecoveryRate(S_FGLQ,S_opt,Q,rho_snr);

end


%% Print success rate of all the tested algorithms
clc;

fprintf("=======OMP_C=======\n")
fprintf("Success Rates: [");
fprintf("%.4f, ", Rate_OMP);
fprintf("]\n");
fprintf("Time: [");
fprintf("%.4f, ", time_OMP);
fprintf("]\n");

fprintf("=======CoSaMP_C=======\n")
fprintf("Success Rates: [");
fprintf("%.4f, ", Rate_CoSaMP);
fprintf("]\n");
fprintf("Time: [");
fprintf("%.4f, ", time_CoSaMP);
fprintf("]\n");

fprintf("=======BP_C=======\n")
fprintf("Success Rates: [");
fprintf("%.4f, ", Rate_BP);
fprintf("]\n");
fprintf("Time: [");
fprintf("%.4f, ", time_BP);
fprintf("]\n");

fprintf("=======Yin_IRLS=======\n")
fprintf("Success Rates: [");
fprintf("%.4f, ", Rate_Yin_IRLS);
fprintf("]\n");
fprintf("Time: [");
fprintf("%.4f, ", time_Yin_IRLS);
fprintf("]\n");

fprintf("=======FISTA=======\n")
fprintf("Success Rates: [");
fprintf("%.4f, ", Rate_FISTA);
fprintf("]\n");
fprintf("Time: [");
fprintf("%.4f, ", time_FISTA);
fprintf("]\n");

fprintf("=======Alg_{GL2}=======\n")
fprintf("Success Rates: [");
fprintf("%.4f, ", Rate_GL2);
fprintf("]\n");
fprintf("Time: [");
fprintf("%.4f, ", time_GL2);
fprintf("]\n");

fprintf("=======Alg_{GL1}=======\n")
fprintf("Success Rates: [");
fprintf("%.4f, ", Rate_GL1);
fprintf("]\n");
fprintf("Time: [");
fprintf("%.4f, ", time_GL1);
fprintf("]\n");

fprintf("=======Alg_{GL1}^F=======\n")
fprintf("Success Rates: [");
fprintf("%.4f, ", Rate_FGL1);
fprintf("]\n");
fprintf("Time: [");
fprintf("%.4f, ", time_FGL1);
fprintf("]\n");

fprintf("=======Alg_{GLQ}=======\n")
fprintf("Success Rates: [");
fprintf("%.4f, ", Rate_GLQ);
fprintf("]\n");
fprintf("Time: [");
fprintf("%.4f, ", time_GLQ);
fprintf("]\n");

fprintf("=======Alg_{GLQ}^F=======\n")
fprintf("Success Rates: [");
fprintf("%.4f, ", Rate_FGLQ);
fprintf("]\n");
fprintf("Time: [");
fprintf("%.4f, ", time_FGLQ);
fprintf("]\n");



