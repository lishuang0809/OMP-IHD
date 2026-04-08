% This code is used to generate Figure 1.
% Coded by Shuang Li, Jan. 2025.

clc; clear; close all;

J = 1000;            % # of samples
N = 64; %64;
L = 128; %128;           % The dictionary is of size N x L
kk = 1:1:N/2;      % sparsity


%alg = [0 1 2 3 4 100 1000 101 102 104 105 106];
rho_snr = 60;   %in dB, used to define a successful recovery.

lenk = length(kk);
Rate_OMP = zeros(1,lenk);
time_OMP = zeros(1,lenk);
Rate_CoSaMP = zeros(1,lenk);
time_CoSaMP = zeros(1,lenk);
Rate_BP = zeros(1,lenk);
time_BP = zeros(1,lenk);
Rate_Yin_IRLS = zeros(1,lenk);
time_Yin_IRLS = zeros(1,lenk);
Rate_FISTA = zeros(1,lenk);
time_FISTA = zeros(1,lenk);
Rate_GL2 = zeros(1,lenk);
time_GL2 = zeros(1,lenk);
Rate_GL1 = zeros(1,lenk);
time_GL1 = zeros(1,lenk);
Rate_FGL1 = zeros(1,lenk);
time_FGL1 = zeros(1,lenk);
Rate_GLQ = zeros(1,lenk);
time_GLQ = zeros(1,lenk);
Rate_FGLQ = zeros(1,lenk);
time_FGLQ = zeros(1,lenk);

for i = 1:lenk
    k = kk(i);
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
clc; close all;

% V_kappa=1:N/2;
t = kk/N;

figure;
fs = 20;
lw = 2;
ms = 10;
num = 8;
nc = 3;

plot(t,Rate_OMP,'k:','linewidth',lw,'MarkerIndices',1:lenk/num:lenk,'markersize', ms,'DisplayName','OMP_C');hold on;
plot(t,Rate_BP,'b-.','linewidth',lw,'MarkerIndices',1:lenk/num:lenk,'markersize', ms,'DisplayName','BP_C');hold on;
plot(t,Rate_CoSaMP,':','Color',[1, 0.5, 0],'linewidth',lw,'MarkerIndices',1:length(t)/num:length(t),'markersize', ms,'DisplayName','CoSaMP_C');hold on;
plot(t,Rate_Yin_IRLS,'g-','linewidth',lw,'MarkerIndices',1:lenk/num:lenk,'markersize', ms,'DisplayName','IRLS_C');hold on;
plot(t,Rate_FISTA,'c--','linewidth',lw,'MarkerIndices',1:lenk/num:lenk,'markersize', ms,'DisplayName','FISTA');hold on;
plot(t,Rate_GL2,'-x','Color',[0.5, 0.2, 1],'linewidth',lw,'MarkerIndices',1:lenk/num:lenk,'markersize', ms,'DisplayName','Alg_{GL2}');hold on;
plot(t,Rate_GL1,'-*','Color',[0.5, 0.8, 1],'linewidth',lw,'MarkerIndices',1:lenk/num:lenk,'markersize', ms,'DisplayName','Alg_{GL1}');hold on;
plot(t,Rate_FGL1,'-o','Color',[0.5, 0, 0.5],'linewidth',lw,'MarkerIndices',1:lenk/num:lenk,'markersize', ms,'DisplayName','Alg_{GL1}^F');hold on;
plot(t,Rate_GLQ,'m-+','linewidth',lw,'MarkerIndices',1:lenk/num:lenk,'markersize', ms,'DisplayName','Alg_{GLQ}');hold on;
plot(t,Rate_FGLQ, 'r-p','linewidth',lw,'MarkerIndices',1:length(t)/num:length(t),'markersize', ms,'DisplayName','Alg_{GLQ}^{F}');hold on;


% plot(t,Rate_CoSaMP,':','Color',[1, 0.5, 0],'linewidth',lw,'MarkerIndices',1:length(t)/num:length(t),'markersize', ms,'DisplayName','CoSaMP_C');hold on;
% plot(t,Rate_FISTA,'c--','linewidth',lw,'MarkerIndices',1:lenk/num:lenk,'markersize', ms,'DisplayName','FISTA');hold on;
% plot(t,Rate_Yin_IRLS,'g-','linewidth',lw,'MarkerIndices',1:lenk/num:lenk,'markersize', ms,'DisplayName','Yin\_IRLS');hold on;
% plot(t,Rate_Alg_GBP,'m-*','linewidth',lw,'MarkerIndices',1:lenk/num:lenk,'markersize', ms,'DisplayName','Alg_{GBP}');hold on;
% plot(t,Rate_super_Fast_Alg_GBP,'r-o','linewidth',lw,'MarkerIndices',1:lenk/num:lenk,'markersize', ms,'DisplayName','Alg_{FGBP}');hold on;
% plot(t,Rate_Alg_l2_l1_1,'-d','Color',[0.5, 0.2, 1],'linewidth',lw,'MarkerIndices',1:length(t)/num:length(t),'markersize', ms,'DisplayName','Alg_{l_2/l_1}^{(1)}');hold on;
% plot(t,Rate_Alg_l2_l1_2, '-p','Color',[0.5, 0, 0.5],'linewidth',lw,'MarkerIndices',1:length(t)/num:length(t),'markersize', ms,'DisplayName','Alg_{l_2/l_1}^{(2)}');hold on;

% set(gca, 'YScale', 'log'); % Logarithmic y-axis
xlabel('$\kappa/N$','Interpreter','LaTex','FontSize',fs)
%ylabel('Rate of successful recovery','Interpreter','LaTex','FontSize',fs)
ylabel('$\rho_{ok}$','Interpreter','LaTex','FontSize',fs)
legend('show','Location', 'southwest','NumColumns', nc,'FontSize',fs-10)
set(gca,'FontSize',fs);
set(gcf, 'Color', 'w');
axis tight

% export_fig 'figs/fig_Data_N64_rho.pdf'; 



figure;
fs = 20;
lw = 2;
ms = 10;
num = 8;
nc = 3;

plot(t,time_OMP,'k:','linewidth',lw,'MarkerIndices',1:lenk/num:lenk,'markersize', ms,'DisplayName','OMP_C');hold on;
plot(t,time_BP,'b-.','linewidth',lw,'MarkerIndices',1:lenk/num:lenk,'markersize', ms,'DisplayName','BP_C');hold on;
plot(t,time_CoSaMP,':','Color',[1, 0.5, 0],'linewidth',lw,'MarkerIndices',1:length(t)/num:length(t),'markersize', ms,'DisplayName','CoSaMP_C');hold on;
plot(t,time_Yin_IRLS,'g-','linewidth',lw,'MarkerIndices',1:lenk/num:lenk,'markersize', ms,'DisplayName','IRLS_C');hold on;
plot(t,time_FISTA,'c--','linewidth',lw,'MarkerIndices',1:lenk/num:lenk,'markersize', ms,'DisplayName','FISTA');hold on;
plot(t,time_GL2,'-x','Color',[0.5, 0.2, 1],'linewidth',lw,'MarkerIndices',1:lenk/num:lenk,'markersize', ms,'DisplayName','Alg_{GL2}');hold on;
plot(t,time_GL1,'-*','Color',[0.5, 0.8, 1],'linewidth',lw,'MarkerIndices',1:lenk/num:lenk,'markersize', ms,'DisplayName','Alg_{GL1}');hold on;
plot(t,time_FGL1,'-o','Color',[0.5, 0, 0.5],'linewidth',lw,'MarkerIndices',1:lenk/num:lenk,'markersize', ms,'DisplayName','Alg_{GL1}^F');hold on;
plot(t,time_GLQ,'m-+','linewidth',lw,'MarkerIndices',1:lenk/num:lenk,'markersize', ms,'DisplayName','Alg_{GLQ}');hold on;
plot(t,time_FGLQ, 'r-p','linewidth',lw,'MarkerIndices',1:length(t)/num:length(t),'markersize', ms,'DisplayName','Alg_{GLQ}^{F}');hold on;


% plot(t,Rate_CoSaMP,':','Color',[1, 0.5, 0],'linewidth',lw,'MarkerIndices',1:length(t)/num:length(t),'markersize', ms,'DisplayName','CoSaMP_C');hold on;
% plot(t,Rate_FISTA,'c--','linewidth',lw,'MarkerIndices',1:lenk/num:lenk,'markersize', ms,'DisplayName','FISTA');hold on;
% plot(t,Rate_Yin_IRLS,'g-','linewidth',lw,'MarkerIndices',1:lenk/num:lenk,'markersize', ms,'DisplayName','Yin\_IRLS');hold on;
% plot(t,Rate_Alg_GBP,'m-*','linewidth',lw,'MarkerIndices',1:lenk/num:lenk,'markersize', ms,'DisplayName','Alg_{GBP}');hold on;
% plot(t,Rate_super_Fast_Alg_GBP,'r-o','linewidth',lw,'MarkerIndices',1:lenk/num:lenk,'markersize', ms,'DisplayName','Alg_{FGBP}');hold on;
% plot(t,Rate_Alg_l2_l1_1,'-d','Color',[0.5, 0.2, 1],'linewidth',lw,'MarkerIndices',1:length(t)/num:length(t),'markersize', ms,'DisplayName','Alg_{l_2/l_1}^{(1)}');hold on;
% plot(t,Rate_Alg_l2_l1_2, '-p','Color',[0.5, 0, 0.5],'linewidth',lw,'MarkerIndices',1:length(t)/num:length(t),'markersize', ms,'DisplayName','Alg_{l_2/l_1}^{(2)}');hold on;

set(gca, 'YScale', 'log'); % Logarithmic y-axis
xlabel('$\kappa/N$','Interpreter','LaTex','FontSize',fs)
ylabel('$T_c$','Interpreter','LaTex','FontSize',fs)
legend('show','Location', 'southeast','NumColumns', nc,'FontSize',fs-10)
set(gca,'FontSize',fs);
set(gcf, 'Color', 'w');
axis tight

% export_fig 'figs/fig_Data_N64_Tc.pdf'; 










