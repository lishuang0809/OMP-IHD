% This code is used to generate phase transition plot in terms of
% sparsity and signal dimension.

% Coded by Shuang Li, February 2024.
% Updated in Jan. 2025.

clc; clear; %close all;

NN = 128:8:256;
J = 100;          % # of samples
kk = 20:5:90;


rho_snr = 60;   % for computing recovery rate

sr_OMP = zeros(length(kk),length(NN));
sr_IRLS = zeros(length(kk),length(NN));
sr_GL2 = zeros(length(kk),length(NN));
sr_GLQ = zeros(length(kk),length(NN));
sr_FGLQ = zeros(length(kk),length(NN));

for in = 1:length(NN)

    for ik = 1:length(kk)
        fprintf('N = %d, kappa = %d\n', NN(in),kk(ik))
        N = NN(in);
        k = kk(ik); % sparsity
        L = 2*N;    % The dictionary is of size N x L
        [X, Q, S_opt] = data_generation(N, L, k, J);  % X = Q*S_opt, ||X(:,j)||_2=1

        % ======= OMP_C ========
        [S_OMP,~] = finalized_Master_4_Lis(X,Q,k,1);
        [sr_OMP(ik,in),~,~] = ComputingRecoveryRate(S_OMP,S_opt,Q,rho_snr);

        % ======= Yin_IRLS ========
        [S_Yin_IRLS,~] = finalized_Master_4_Lis(X,Q,k,4);
        [sr_IRLS(ik,in),~,~] = ComputingRecoveryRate(S_Yin_IRLS,S_opt,Q,rho_snr);

        % ======= Alg_{GL2} ========
        [S_GL2,~] = finalized_Master_4_Lis(X,Q,k,6);
        [sr_GL2(ik,in),~,~] = ComputingRecoveryRate(S_GL2,S_opt,Q,rho_snr);

        % ======= Alg_{GLQ} ========
        [S_GLQ,~] = finalized_Master_4_Lis(X,Q,k,9);
        [sr_GLQ(ik,in),~,~] = ComputingRecoveryRate(S_GLQ,S_opt,Q,rho_snr);

        % ======= Alg_{GLQ}^F ========
        [S_FGLQ,~] = finalized_Master_4_Lis(X,Q,k,10);
        [sr_FGLQ(ik,in),~,~] = ComputingRecoveryRate(S_FGLQ,S_opt,Q,rho_snr);


    end
    fprintf('Progress: %d/%d (%.2f%%)\n', in, length(NN), (in / length(NN)) * 100);

end


%% Show results

figure(1);clf;
fs=20;
imagesc(NN,kk,sr_OMP);axis xy
xlabel('$N$','Interpreter','LaTex','FontSize',fs)
ylabel('Sparsity ($\kappa$)','Interpreter','LaTex','FontSize',fs);
title('Rate of successful recovery (OMP_C)','FontSize',fs);
set(gcf, 'Color', 'w');
colormap pink;
colorbar; clim([0 1]);
set(gca,'FontSize',fs);
%export_fig 'figs/test_phasetransition_OMP_C.pdf';

figure(2);clf;
fs=20;
imagesc(NN,kk,sr_IRLS);axis xy
xlabel('$N$','Interpreter','LaTex','FontSize',fs)
ylabel('Sparsity ($\kappa$)','Interpreter','LaTex','FontSize',fs);
title('Rate of successful recovery (IRLS_C)','FontSize',fs);
set(gcf, 'Color', 'w');
colormap pink;
colorbar; clim([0 1]);
set(gca,'FontSize',fs);
%export_fig 'figs/test_phasetransition_IRLS_C.pdf';

figure(3);clf;
fs=20;
imagesc(NN,kk,sr_GL2);axis xy
xlabel('$N$','Interpreter','LaTex','FontSize',fs)
ylabel('Sparsity ($\kappa$)','Interpreter','LaTex','FontSize',fs);
title('Rate of successful recovery (Alg_{GL2})','FontSize',fs);
set(gcf, 'Color', 'w');
colormap pink;
colorbar; clim([0 1]);
set(gca,'FontSize',fs);
%export_fig 'figs/test_phasetransition_GL2.pdf';

figure(4);clf;
fs=20;
imagesc(NN,kk,sr_GLQ);axis xy
xlabel('$N$','Interpreter','LaTex','FontSize',fs)
ylabel('Sparsity ($\kappa$)','Interpreter','LaTex','FontSize',fs);
title('Rate of successful recovery (Alg_{GLQ})','FontSize',fs);
set(gcf, 'Color', 'w');
colormap pink;
colorbar; clim([0 1]);
set(gca,'FontSize',fs);
%export_fig 'figs/test_phasetransition_GLQ.pdf';


figure(5);clf;
fs=20;
imagesc(NN,kk,sr_FGLQ);axis xy
xlabel('$N$','Interpreter','LaTex','FontSize',fs)
ylabel('Sparsity ($\kappa$)','Interpreter','LaTex','FontSize',fs);
title('Rate of successful recovery (Alg_{GLQ}^F)','FontSize',fs);
set(gcf, 'Color', 'w');
colormap pink;
colorbar; clim([0 1]);
set(gca,'FontSize',fs);
%export_fig 'figs/test_phasetransition_FGLQ.pdf';


figure(6);clf;
fs=20;
imagesc(NN,kk,sr_FGLQ - sr_IRLS);axis xy
xlabel('$N$','Interpreter','LaTex','FontSize',fs)
ylabel('Sparsity ($\kappa$)','Interpreter','LaTex','FontSize',fs);
title('Alg_{GLQ}^F - IRLS_C','FontSize',fs);
set(gcf, 'Color', 'w');
colormap pink;
colorbar; %clim([0 1]);
set(gca,'FontSize',fs);
%export_fig 'figs/test_phasetransition_diff.pdf';









