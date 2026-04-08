% This is to plot (R,T) v.s sigma_e
                
function[V_R,V_T] = Run_finalized_Master_4_Lis(Q,S_opt,kappa,sigma_e,alg)
        % The input data can be generated with
        %    [Q,S_opt]=generating_data(N,L,J,kappa);
        %        [Q,S_opt,Q_condi]=gen_data_condQ(N,L,J,kappa,condi)
        
  [N,L]= size(Q);                [~,J] = size(S_opt);
  % sigma_e = 0;      
       % alg = 104;   % [1 2 3 4 5 6 7 8 9 10]; 
                                               % Algorithms to be examined     
  % Setting for computing ratio of successful sparse recovery         
  
           rho_snr = 60;   %in dB                                             
   
 % Generating the clean signals and the noises
   X = Q*S_opt;          E = diag(X'*X);      E = sqrt(E);    % J x 1
   S_opt = S_opt*inv(diag(E));  
   X_opt = Q*S_opt;        % Such S_opt makes  ||X(:,j)||_2 = 1, for all j
   
   E = randn(N,J);   E = normc(E);   
   E = E*sigma_e;    % Note ||E(:,j)||_2 = sigma_e
   
   X = X_opt + E;        % The signal-to-noise ratio for each j is  1/sigma_e
                     % rho_snr = 20*log_10(1/sigma_e)   dB
 %====================================================================%   
  V_R=[]; V_T=[];            %V_snr=[];  V_snr_ref = [];
  
 for k=1:length(alg)
                     [S,T] = finalized_Master_4_Lis(X,Q,kappa,alg(k));    
% The original S, Z                   
   [R,~] = ComputingRecoveryRate(S,S_opt,Q,rho_snr)  ; 
                 % snr - 1 x J
  V_R = [V_R R];           V_T = [V_T T];        % 1 x length(alg) 
  %V_snr = [V_snr sum(snr)/J];  % averaged snr over the J samples
 
 end
 % M_R = [M_R;V_R];   M_T = [M_T; V_T];   M_snr = [M_snr; V_snr];
  
   % M_R - length(V_kappa) x length(alg)


 %save data_4_Lis V_R V_T
 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% F1=figure;
% set(F1,'PaperUnits', 'centimeters'); 
% set(F1,'Position',[0 0 600 400]);
% 
% 
% F11=subplot(2,1,1);
% set(F11,'LineWidth',1,'FontName','Times New Roman','FontSize',12,'FontWeight','Bold');
% hold on;
% 
% plot(kappa,V_4(:,1), '-v')   %  虚线
% hold on;
% plot(kappa,V_7(:,1), '-*')   % 点线
% hold on;
% plot(kappa,V_8(:,1), '-x')   % 点划�?
% 
% axis([1 32  0 1]) 
%  xlabel('(a)')
% % ylabel('x(t),  x_1(t)');
% 
% 
% F12=subplot(2,1,2);
% set(F12,'LineWidth',1,'FontName','Times New Roman','FontSize',12,'FontWeight','Bold');
% hold on;
% 
% plot(kappa,V_4(:,2), '-v')   %  虚线
% hold on;
% plot(kappa,V_7(:,2), '-*')   % 点线
% hold on;
% plot(kappa,V_8(:,2), '-x')   % 点划�?
% 
% axis([1 32  0 20]) 
%  xlabel('(b)')
% % ylabel('x(t),  x_1(t)');
% 
% %print('fig-K_varying','-deps')   %After the figure shown, it can be save as a
% 
% saveas(gcf,"fig-K_varying")            % saved in fig-format
% saveas(gcf,"fig-K_varying",'epsc')     % saved in eps-format
% 

%%%%%%