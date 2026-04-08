clc; clear;
N = 64;
L = 128;     % The dictionary is of size N x L
k = 20;          % sparsity
J = 10;      % 1000
[X, Q, S_opt] = data_generation(N, L, k, J);  % X = Q*S_opt
norm(X-Q*S_opt,'fro')

rho_snr = 60;   %in dB

% ======= OMP_C ========
[S_OMP,time_OMP] = finalized_Master_4_Lis(X,Q,k,1);
[Rate_OMP,~,~] = ComputingRecoveryRate(S_OMP,S_opt,Q,rho_snr);

% ======= BP_C ========
[S_BP,time_BP] = finalized_Master_4_Lis(X,Q,k,3);
[Rate_BP,~,~] = ComputingRecoveryRate(S_BP,S_opt,Q,rho_snr);

% ======= Alg_GL2 ========
[S_Alg_GL2,time_Alg_GL2] = finalized_Master_4_Lis(X,Q,k,6);
[Rate_Alg_GL2,~,~] = ComputingRecoveryRate(S_Alg_GL2,S_opt,Q,rho_snr);

%=======Alg_GL1=======
[S_Alg_GL1,time_Alg_GL1] = finalized_Master_4_Lis(X,Q,k,7);
[Rate_Alg_GL1,~,~] = ComputingRecoveryRate(S_Alg_GL1,S_opt,Q,rho_snr);


clc;
fprintf("=======OMP=======\n")
fprintf("Success Rate: %.4f\n",Rate_OMP)
fprintf("Running time: %.4f\n",time_OMP)

fprintf("=======BP=======\n")
fprintf("Success Rate: %.4f\n",Rate_BP)
fprintf("Running time: %.4f\n",time_BP)

fprintf("=======Alg_GL2=======\n")
fprintf("Success Rate: %.4f\n",Rate_Alg_GL2)
fprintf("Running time: %.4f\n",time_Alg_GL2)

fprintf("=======Alg_GL1=======\n")
fprintf("Success Rate: %.4f\n",Rate_Alg_GL1)
fprintf("Running time: %.4f\n",time_Alg_GL1)


