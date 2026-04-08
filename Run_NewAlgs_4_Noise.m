% The codes below is to deal with x = Qs + e, where e = Gamma*alpha with 
%      Gamma^T*Gamma = I_N_e  and N_e <= N, based on our proposed
%      characterization: s = T_0 (x-e) + W*z 
%                          = T_0*x + [W -T_0*Gamma]*[z^T alpha^T]^T
                
function[V_R,V_T,V_RE] = Run_NewAlgs_4_Noise(Q,S_opt,kappa,N_e,sigma_e,V_alg)
        % The input data can be generated with
        %    [Q,S_opt]=generating_data(N,L,J,kappa);
        %    V_alp = [104 107 108 109]
        
  [N,L]= size(Q);                [~,J] = size(S_opt);
     
  rho_snr = 60;   %in dB                                             
   
 % Generating the clean signals and the noises
   X = Q*S_opt;          E = diag(X'*X);      E = sqrt(E);    % J x 1
   S_opt = S_opt*inv(diag(E));  
   X_opt = Q*S_opt;        % Such S_opt makes  ||X(:,j)||_2 = 1, for all j
   
   Gamma = randn(N);   [Gamma,~,~]=svd(Gamma); Gamma = Gamma(:,1:N_e);
   V_alpha = randn(N_e,J); V_alpha = normc(V_alpha); 
   V_alpha = V_alpha*sigma_e;    
   E = Gamma*V_alpha;          % Note ||E(:,j)||_2 = sigma_e
   
   X = X_opt + E;        % The signal-to-noise ratio for each j is  1/sigma_e
                         % rho_snr = 20*log_10(1/sigma_e)   dB
          %============================================%   
 % Preparation 
  [U,Sig,V]=svd(Q);  Sig=Sig(:,1:N);        % Q = U[Sig 0]V^T  
  V_1 = V(:,1:N);                % Q_bar = V_1*V_1';
  T_0 = V(:,1:N)/Sig*U';         %SL: was V(:,1:N)*inv(Sig)*U';   
  W = V(:,N+1:L);     
  W_c = Sig\U';                  %SL: was inv(Sig)*U';         
        %T_Q = T_0*Q;
  S_0 = T_0*X;
  
  if sigma_e==0
          W_bar = W;
  else

      
          W_bar = [W -T_0*Gamma];
  end
         %==================================%
 
    %====================================================================%   
  V_R=[]; V_T=[]; V_RE = [];           %V_snr=[];  V_snr_ref = [];
  
 for k=1:length(V_alg)
                 % [S,T] = finalized_Master_4_Lis(X,Q,kappa,V_alg(k));  
        [S,T] = NewAlgs_4_Noise(Q,S_opt,X,kappa,W_bar,S_0,V_1,V_alg(k),Gamma,W,T_0);
% The original S, Z                   
   [R,~] = ComputingRecoveryRate(S,S_opt,Q,rho_snr)  ; 
   RE = norm(S - S_opt,'fro')/norm(S_opt,'fro'); % compute relative error.
                 % snr - 1 x J
  V_R = [V_R R];       V_T = [V_T T];      V_RE = [V_RE RE];       % 1 x length(alg) 
  %V_snr = [V_snr sum(snr)/J];  % averaged snr over the J samples
 
 end

 
end      
         
 