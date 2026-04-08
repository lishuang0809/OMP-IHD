% The codes below is to deal with x = Qs + e, where e = Gamma*alpha with 
%      Gamma^T*Gamma = I_N_e  and N_e <= N, based on our proposed
%      characterization: s = T_0 (x-e) + W*z 
%                          = T_0*x + [W - T_0*Gamma]*[z^T alpha^T]^T
                
function[S_kappa,T_c,X_est] = NewAlgs_4_Noise_updated(Q,S_opt,X,kappa,W_bar,S_0,V_1,alg)
 
 [~,J] = size(S_opt);        S_kappa=S_opt*0;    % estimated sparse vectors
 X_est = X*0;  % estimation of X_opt= Q*S_opt
 %==============================================
 
 if alg==100         % IRLS_alt.m  - NOT in (iterative) greedy manner !!!
    tic,                       
     for j=1:J
                      q=1/2; 
        [s, ~] = IRLS_alt3(S_0(:,j),W_bar,q)  ;
        %===================================================
        s_kappa = H_kappaOnly(s,kappa)  ; I_c = find(abs(s_kappa)==0); 
        D = W_bar(I_c,:);    z_bar = -inv(D'*D)*D'*S_0(I_c,j);
        s = S_0(:,j)+W_bar*z_bar;
        S_kappa(:,j) = s; % estimate of the kappa-sparse vector by Our Oracle
                % estimator
        X_est(:,j) = Q*s;  
                  % esti. error of the classical Oracle estimator
        %=============================================================          
     end
  T_c = toc;   % the time taken in seconds since tic      
end 
%===========================================     

 
if alg==104     % 
  tic,                       
     for j=1:J
                      q=1/2; 
        s = Yin_IRLS_cla3(X(:,j),Q,q)  ;
        %===================================================
        s_kappa = H_kappaOnly(s,kappa)  ; I = find(abs(s_kappa)>0); 
        D = Q(:,I);    beta = inv(D'*D)*D'*X(:,j);
        S_kappa(I,j) =beta;     % estimate of the kappa-sparse vector
        X_est(:,j) = D*beta;  
                  % esti. error of the classical Oracle estimator
        %=============================================================  
     end
  T_c = toc;   % the time taken in seconds since tic      
end 
%===========================================        
 %==============================================
if alg==107     %  the original Alg_GBP - very slow one! 
 tic,   
 for j=1:J
   s = Alg_GBP(W_bar,V_1,S_0(:,j),kappa);
        %===================================================
        s_kappa = H_kappaOnly(s,kappa)  ; I_c = find(abs(s_kappa)==0); 
        D = W_bar(I_c,:);    z_bar = -inv(D'*D)*D'*S_0(I_c,j);
        s = S_0(:,j)+W_bar*z_bar;
        S_kappa(:,j) = s; % estimate of the kappa-sparse vector by Our Oracle
                % estimator
        X_est(:,j) = Q*s;  
                  % esti. error of the classical Oracle estimator
        %============================================================= 
 end
 T_c = toc;   % the time taken in seconds since tic   
end
%============================================        


 %==============================================
if alg==108     %  the original Alg_GBP - very slow one! 
   N_iter=5; 
 tic,   
 for j=1:J
   % s = Alg_GBP_rev(W,V_1,S_0(:,j),kappa);
   s = Alg_FGBP(W_bar,V_1,S_0(:,j),kappa,N_iter);
        %===================================================
        s_kappa = H_kappaOnly(s,kappa)  ; I_c = find(abs(s_kappa)==0); 
        D = W_bar(I_c,:);    z_bar = -inv(D'*D)*D'*S_0(I_c,j);
        s = S_0(:,j)+W_bar*z_bar;
        S_kappa(:,j) = s; % estimate of the kappa-sparse vector by Our Oracle
                % estimator
        X_est(:,j) = Q*s;  
                  % esti. error of the classical Oracle estimator
        %============================================================= 
 end
 T_c = toc;   % the time taken in seconds since tic   
end
%============================================
 
%=================================================
if alg==109     % psi(t) = [|t|^2+epsilon]^{q/2} - l_q-norm

    tic,                       
    for j=1:J
       s = Alg_GLQ(W_bar,V_1,S_0(:,j),kappa);
        %===================================================
        s_kappa = H_kappaOnly(s,kappa)  ; I_c = find(abs(s_kappa)==0); 
        D = W_bar(I_c,:);    z_bar = -inv(D'*D)*D'*S_0(I_c,j);
        s = S_0(:,j)+W_bar*z_bar;
        S_kappa(:,j) = s; % estimate of the kappa-sparse vector by Our Oracle
                % estimator
        X_est(:,j) = Q*s;  
                  % esti. error of the classical Oracle estimator
        %============================================================= 
    end
    T_c=toc;   
end
%===================================================================

%=================================================
if alg==110     % psi(t) = |t|/[|t} + epsil]
          N_iter=2; 
    tic,                       
    for j=1:J
        %s=Alg_FGDM2(W,V_1,S_0(:,j),kappa,N_iter);
        s = Alg_FGLQ(W_bar,V_1,S_0(:,j),kappa,N_iter);
        %===================================================
        s_kappa = H_kappaOnly(s,kappa)  ; I_c = find(abs(s_kappa)==0); 
        D = W_bar(I_c,:);    z_bar = -inv(D'*D)*D'*S_0(I_c,j);
        s = S_0(:,j)+W_bar*z_bar;
        S_kappa(:,j) = s; % estimate of the kappa-sparse vector by Our Oracle
                % estimator
        X_est(:,j) = Q*s;  
                  % esti. error of the classical Oracle estimator
        %============================================================= 
    end
    T_c=toc;   
end
%===================================================================

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