clc; clear; rng(50); 

N = 128;
J = 100; %100;         % # of samples
L = 256;
kappa = 30;         % sparsity

[Q,S_opt]=generating_data(N,L,J,kappa);
%V_alg = [104 107 108 109 110]; % algs: IRLS, GL1, FGL1, GLQ, FGLQ
%V_alg = [104 109 110 111]; % algs: IRLS, GLQ, FGLQ, Proj_FGLQ,  
V_alg = [104 109 110 ]; % algs: IRLS, GLQ, FGLQ 
lena = length(V_alg);

N_e = 20; 
Sigma_e = 0.00:0.01:0.1;    %0.0125;
lens = length(Sigma_e);
Rate = zeros(lens,lena);
Time = zeros(lens,lena);
RError = zeros(lens,lena);

for i = 1:lens
    sigma_e = Sigma_e(i);
    fprintf('Sigma = %.3f\n', sigma_e)
    [R,T_c,RE] = Run_NewAlgs_4_Noise(Q,S_opt,kappa,N_e,sigma_e,V_alg);
    Rate(i,:) = R;
    Time(i,:) = T_c;
    RError(i,:) = RE;
end

%%
Rate_Yin_IRLS = Rate(:,1);
Rate_GLQ = Rate(:,2);
Rate_FGLQ = Rate(:,3);
% Rate_GLQ_proj = Rate(:,4);
%Rate_RPCA = Rate(:,5);
time_Yin_IRLS = Time(:,1);
time_GLQ = Time(:,2);
time_FGLQ = Time(:,3);
% time_GLQ_proj = Time(:,4);
%time_RPCA = Time(:,5);
error_Yin_IRLS = RError(:,1);
error_GLQ = RError(:,2);
error_FGLQ = RError(:,3);
% error_GLQ_proj = RError(:,4);
%error_RPCA = RError(:,5);
%% Show results
figure(1); clf;
fs = 20;
lw = 2;
ms = 10;
num = lens;
nc = 3;
%plot(Sigma_e,Rate_OMP,'k:','linewidth',lw,'MarkerIndices',1:lens/num:lens,'markersize', ms,'DisplayName','OMP_C');hold on;
plot(Sigma_e,Rate_Yin_IRLS,'g-o','linewidth',lw,'MarkerIndices',1:lens/num:lens,'markersize', ms,'DisplayName','IRLS_C');hold on;
%plot(Sigma_e,Rate_FISTA,'c--','linewidth',lw,'MarkerIndices',1:lens/num:lens,'markersize', ms,'DisplayName','FISTA');hold on;
%plot(Sigma_e,Rate_GL2,'-x','Color',[0.5, 0.2, 1],'linewidth',lw,'MarkerIndices',1:lens/num:lens,'markersize', ms,'DisplayName','Alg_{GL2}');hold on;
plot(Sigma_e,Rate_GLQ,'m-+','linewidth',lw,'MarkerIndices',1:lens/num:lens,'markersize', ms,'DisplayName','Alg_{GLQ}');hold on;
plot(Sigma_e,Rate_FGLQ, 'r-p','linewidth',lw,'MarkerIndices',1:lens/num:lens,'markersize', ms,'DisplayName','Alg_{GLQ}^{F}');hold on;
% plot(Sigma_e,Rate_GLQ_proj,'b-x','linewidth',lw,'MarkerIndices',1:lens/num:lens,'markersize', ms,'DisplayName','Alg^F_{GLQ-Proj}');hold on;
%plot(Sigma_e,Rate_RPCA,'k-*','linewidth',lw,'MarkerIndices',1:lens/num:lens,'markersize', ms,'DisplayName','RPCA');hold on;



xlabel('$\sigma$','Interpreter','LaTex','FontSize',fs)
ylabel('$\rho_{ok}$','Interpreter','LaTex','FontSize',fs)
%ylabel('$\|X - \hat{X}\|_F/\|X\|_F$','Interpreter','LaTex','FontSize',fs)
legend('show','Location', 'northeast','NumColumns', nc,'FontSize',fs-8)  %northeast
%set(gca, 'YScale', 'log');
set(gca,'FontSize',fs);
set(gcf, 'Color', 'w');

%export_fig 'figs/test_noise_rho_N128_k30.pdf';


% figure(2); clf;
% fs = 20;
% lw = 2;
% ms = 10;
% num = lens;
% nc = 3;
% %plot(Sigma_e,Rate_OMP,'k:','linewidth',lw,'MarkerIndices',1:lens/num:lens,'markersize', ms,'DisplayName','OMP_C');hold on;
% plot(Sigma_e,error_Yin_IRLS,'g-o','linewidth',lw,'MarkerIndices',1:lens/num:lens,'markersize', ms,'DisplayName','IRLS_C');hold on;
% %plot(Sigma_e,Rate_FISTA,'c--','linewidth',lw,'MarkerIndices',1:lens/num:lens,'markersize', ms,'DisplayName','FISTA');hold on;
% %plot(Sigma_e,Rate_GL2,'-x','Color',[0.5, 0.2, 1],'linewidth',lw,'MarkerIndices',1:lens/num:lens,'markersize', ms,'DisplayName','Alg_{GL2}');hold on;
% plot(Sigma_e,error_GLQ,'m-+','linewidth',lw,'MarkerIndices',1:lens/num:lens,'markersize', ms,'DisplayName','Alg_{GLQ}');hold on;
% plot(Sigma_e,error_FGLQ, 'r-p','linewidth',lw,'MarkerIndices',1:lens/num:lens,'markersize', ms,'DisplayName','Alg_{GLQ}^{F}');hold on;
% % plot(Sigma_e,error_GLQ_proj,'b-x','linewidth',lw,'MarkerIndices',1:lens/num:lens,'markersize', ms,'DisplayName','Alg^F_{GLQ-Proj}');hold on;
% %plot(Sigma_e,error_RPCA,'k-*','linewidth',lw,'MarkerIndices',1:lens/num:lens,'markersize', ms,'DisplayName','RPCA');hold on;
% 
% 
% 
% xlabel('$\sigma$','Interpreter','LaTex','FontSize',fs)
% %ylabel('$\rho_{ok}$','Interpreter','LaTex','FontSize',fs)
% ylabel('$\|S - \hat{S}\|_F/\|S\|_F$','Interpreter','LaTex','FontSize',fs)
% legend('show','Location', 'northeast','NumColumns', nc,'FontSize',fs-8)  %northeast
% set(gca, 'YScale', 'log');
% set(gca,'FontSize',fs);
% set(gcf, 'Color', 'w');
% 
% %export_fig 'figs/test_noise_rho_N128_k80.pdf';



figure(2); clf;
fs = 20;
lw = 2;
ms = 10;
num = lens;
nc = 3;
%plot(Sigma_e,time_OMP,'k:','linewidth',lw,'MarkerIndices',1:lens/num:lens,'markersize', ms,'DisplayName','OMP_C');hold on;
plot(Sigma_e,time_Yin_IRLS,'g-o','linewidth',lw,'MarkerIndices',1:lens/num:lens,'markersize', ms,'DisplayName','IRLS_C');hold on;
%plot(Sigma_e,time_FISTA,'c--','linewidth',lw,'MarkerIndices',1:lens/num:lens,'markersize', ms,'DisplayName','FISTA');hold on;
%plot(Sigma_e,time_GL2,'-x','Color',[0.5, 0.2, 1],'linewidth',lw,'MarkerIndices',1:lens/num:lens,'markersize', ms,'DisplayName','Alg_{GL2}');hold on;
plot(Sigma_e,time_GLQ,'m-+','linewidth',lw,'MarkerIndices',1:lens/num:lens,'markersize', ms,'DisplayName','Alg_{GLQ}');hold on;
plot(Sigma_e,time_FGLQ, 'r-p','linewidth',lw,'MarkerIndices',1:lens/num:lens,'markersize', ms,'DisplayName','Alg_{GLQ}^{F}');hold on;
% plot(Sigma_e,time_GLQ_proj,'b-x','linewidth',lw,'MarkerIndices',1:lens/num:lens,'markersize', ms,'DisplayName','Alg^F_{GLQ-Proj}');hold on;
%plot(Sigma_e,time_RPCA,'k-*','linewidth',lw,'MarkerIndices',1:lens/num:lens,'markersize', ms,'DisplayName','RPCA');hold on;



xlabel('$\sigma$','Interpreter','LaTex','FontSize',fs)
ylabel('$T_c$','Interpreter','LaTex','FontSize',fs)
legend('show','Location', 'northwest','NumColumns', nc,'FontSize',fs-8)  %northeast
%set(gca, 'YScale', 'log');
set(gca,'FontSize',fs);
set(gcf, 'Color', 'w');

%export_fig 'figs/test_noise_Tc_N128_k30.pdf';
