
% This code is used to compare the proposed algorithms with xxx algorithms
% in terms of Q with different condition numbers.

% Coded by Shuang Li, Jan. 2025.


clc; clear; %close all;

N = 256;
J = 10;         % # of samples
L = 512;
k = 25;   %round(0.2*N);          % sparsity

rho_snr = 60;   % for computing recovery rate

CC_Q = [5 10 50 100 500 1e3 5e3 1e4 5e4 1e5];
lens = length(CC_Q);
Rate_OMP = zeros(1,lens);
time_OMP = zeros(1,lens);
% Rate_CoSaMP = zeros(1,lenk);
% time_CoSaMP = zeros(1,lenk);
% Rate_BP = zeros(1,lenk);
% time_BP = zeros(1,lenk);
Rate_Yin_IRLS = zeros(1,lens);
time_Yin_IRLS = zeros(1,lens);
% Rate_FISTA = zeros(1,lens);
% time_FISTA = zeros(1,lens);
Rate_GL2 = zeros(1,lens);
time_GL2 = zeros(1,lens);
% Rate_GL1 = zeros(1,lenk);
% time_GL1 = zeros(1,lenk);
% Rate_FGL1 = zeros(1,lenk);
% time_FGL1 = zeros(1,lenk);
Rate_GLQ = zeros(1,lens);
time_GLQ = zeros(1,lens);
Rate_FGLQ = zeros(1,lens);
time_FGLQ = zeros(1,lens);


for i = 1:lens
    C_Q = CC_Q(i);
    fprintf('C_Q = %.3f\n', CC_Q(i))

    %[Xclean, Q, S_opt] = data_generation(N, L, k, J);  % X = Q*S_opt, ||X(:,j)||_2=1

    [Q, S_opt, Q_condi] = gen_data_condQ(N, L, J, k, C_Q);   
    X = Q_condi*S_opt;  % cond(Q_condi) = condi, ||Q_condi||^2_F=L
    E = diag(X'*X);      E = sqrt(E);
    S_opt = S_opt*inv(diag(E));  
    X = Q_condi*S_opt;        % Such S_opt makes  ||X(:,j)||_2 = 1, for all j

    % ======= OMP_C ========
    fprintf('Run OMP\n')
    [S_OMP,time_OMP(i)] = finalized_Master_4_Lis(X,Q_condi,k,1);
    [Rate_OMP(i),~,~] = ComputingRecoveryRate(S_OMP,S_opt,Q_condi,rho_snr);

    % ======= CoSaMP_C ========
    % [S_CoSaMP,time_CoSaMP(i)] = finalized_Master_4_Lis(X,Q,k,2);
    % [Rate_CoSaMP(i),~,~] = ComputingRecoveryRate(S_CoSaMP,S_opt,Q,rho_snr);

    % ======= BP_C ========
    % [S_BP,time_BP(i)] = finalized_Master_4_Lis(X,Q,k,3);
    % [Rate_BP(i),~,~] = ComputingRecoveryRate(S_BP,S_opt,Q,rho_snr);

    % ======= Yin_IRLS ========
    fprintf('Run Yin_IRLS\n')
    [S_Yin_IRLS,time_Yin_IRLS(i)] = finalized_Master_4_Lis(X,Q_condi,k,4);
    [Rate_Yin_IRLS(i),~,~] = ComputingRecoveryRate(S_Yin_IRLS,S_opt,Q_condi,rho_snr);

    % ======= FISTA ========
    % [S_FISTA,time_FISTA(i)] = finalized_Master_4_Lis(X,Q,k,5);
    % [Rate_FISTA(i),~,~] = ComputingRecoveryRate(S_FISTA,S_opt,Q,rho_snr);


    % ======= Alg_{GL2} ========
    fprintf('Run Alg_{GL2}\n')
    [S_GL2,time_GL2(i)] = finalized_Master_4_Lis(X,Q_condi,k,6);
    [Rate_GL2(i),~,~] = ComputingRecoveryRate(S_GL2,S_opt,Q_condi,rho_snr);

    %======= Alg_{GL1} =======
    % [S_GL1,time_GL1(i)] = finalized_Master_4_Lis(X,Q,k,7);
    % [Rate_GL1(i),~,~] = ComputingRecoveryRate(S_GL1,S_opt,Q,rho_snr);

    %======= Alg_{GL1}^F =======
    % [S_FGL1,time_FGL1(i)] = finalized_Master_4_Lis(X,Q,k,8);
    % [Rate_FGL1(i),~,~] = ComputingRecoveryRate(S_FGL1,S_opt,Q,rho_snr);

    % ======= Alg_{GLQ} ========
    fprintf('Run Alg_{GLQ}\n')
    [S_GLQ,time_GLQ(i)] = finalized_Master_4_Lis(X,Q_condi,k,9);
    [Rate_GLQ(i),~,~] = ComputingRecoveryRate(S_GLQ,S_opt,Q_condi,rho_snr);

    % ======= Alg_{GLQ}^F ========
    fprintf('Run Alg_{GLQ}^F\n')
    [S_FGLQ,time_FGLQ(i)] = finalized_Master_4_Lis(X,Q_condi,k,10);
    [Rate_FGLQ(i),~,~] = ComputingRecoveryRate(S_FGLQ,S_opt,Q_condi,rho_snr);


    fprintf('Progress: %d/%d (%.2f%%)\n', i, lens, (i / lens) * 100);
end

%% Show results
figure(1); clf;
fs = 20;
lw = 2;
ms = 10;
num = lens;
nc = 3;
plot(CC_Q,Rate_OMP,'k:','linewidth',lw,'MarkerIndices',1:lens/num:lens,'markersize', ms,'DisplayName','OMP_C');hold on;
plot(CC_Q,Rate_Yin_IRLS,'g-','linewidth',lw,'MarkerIndices',1:lens/num:lens,'markersize', ms,'DisplayName','IRLS_C');hold on;
% plot(CC_Q,Rate_FISTA,'c--','linewidth',lw,'MarkerIndices',1:lens/num:lens,'markersize', ms,'DisplayName','FISTA');hold on;
plot(CC_Q,Rate_GL2,'-x','Color',[0.5, 0.2, 1],'linewidth',lw,'MarkerIndices',1:lens/num:lens,'markersize', ms,'DisplayName','Alg_{GL2}');hold on;
plot(CC_Q,Rate_GLQ,'m-+','linewidth',lw,'MarkerIndices',1:lens/num:lens,'markersize', ms,'DisplayName','Alg_{GLQ}');hold on;
plot(CC_Q,Rate_FGLQ, 'r-p','linewidth',lw,'MarkerIndices',1:lens/num:lens,'markersize', ms,'DisplayName','Alg_{GLQ}^{F}');hold on;

sq = 1:5;
CC_Qshort = 10.^sq;
%xticks(CC_Qshort);   
xticks(CC_Q);  
%xticklabels(arrayfun(@num2str, CC_Q, 'UniformOutput', false)); 


xlabel('$\mathcal{C}_Q$','Interpreter','LaTex','FontSize',fs)
ylabel('$\rho_{ok}$','Interpreter','LaTex','FontSize',fs)
legend('show','Location', 'northeast','NumColumns', nc,'FontSize',fs-8)  %northeast
set(gca, 'XScale', 'log');
set(gca,'FontSize',fs);
set(gcf, 'Color', 'w');

%export_fig 'figs/test_condQ_rho_N512.pdf';


figure(2); clf;
fs = 20;
lw = 2;
ms = 10;
num = lens;
nc = 3;
plot(CC_Q,time_OMP,'k:','linewidth',lw,'MarkerIndices',1:lens/num:lens,'markersize', ms,'DisplayName','OMP_C');hold on;
plot(CC_Q,time_Yin_IRLS,'g-','linewidth',lw,'MarkerIndices',1:lens/num:lens,'markersize', ms,'DisplayName','IRLS_C');hold on;
%plot(CC_Q,time_FISTA,'c--','linewidth',lw,'MarkerIndices',1:lens/num:lens,'markersize', ms,'DisplayName','FISTA');hold on;
plot(CC_Q,time_GL2,'-x','Color',[0.5, 0.2, 1],'linewidth',lw,'MarkerIndices',1:lens/num:lens,'markersize', ms,'DisplayName','Alg_{GL2}');hold on;
plot(CC_Q,time_GLQ,'m-+','linewidth',lw,'MarkerIndices',1:lens/num:lens,'markersize', ms,'DisplayName','Alg_{GLQ}');hold on;
plot(CC_Q,time_FGLQ, 'r-p','linewidth',lw,'MarkerIndices',1:lens/num:lens,'markersize', ms,'DisplayName','Alg_{GLQ}^{F}');hold on;

%xticks(CC_Qshort);   
xticks(CC_Q);   
%xticklabels(arrayfun(@num2str, CC_Q, 'UniformOutput', false));


xlabel('$\mathcal{C}_Q$','Interpreter','LaTex','FontSize',fs)
ylabel('$T_c$','Interpreter','LaTex','FontSize',fs)
legend('show','Location', 'northwest','NumColumns', nc,'FontSize',fs-10)  %northeast
set(gca, 'XScale', 'log');
% set(gca, 'YScale', 'log');
set(gca,'FontSize',fs);
set(gcf, 'Color', 'w');


%export_fig 'figs/test_condQ_Tc_N2048.pdf';
