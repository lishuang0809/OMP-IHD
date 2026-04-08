% This prog. is to  
%           i) find the kappa-sparse vector     s_hat = H_kappa(s)
%                        with the kappa largest  entriies in mangnitude  s
%          ii) check if I_s_hat = I_s_opt
                
function[RecoverRate,V_snr,Err_index] = ComputingRecoveryRate(S,S_opt,Q,rho_snr)          
            % S_opt   - the ideal kappa-sparse 
            % S       - vectors to be checked
            % kappa   - sparssity level given  - NOT needed below
% Initializing
                      [~,J]=size(S);  

  Err_index = zeros(1,J);  % Err_index(j)=1, S_kappa(:,j) - successful !
  Err_LAMBDA = zeros(1,J);  V_snr=zeros(1,J);
for j=1:J
   % The support of s_opt 
     I_opt = find(abs(S_opt(:,j))>0); % suuport of s_opt 
     L_opt = length(I_opt);           % L_opt <= kappa
    
   % The corresponding support of s to be examined
     [~,Ind]=sort(-abs(S(:,j)));   % sorting the entries 
     Ind_p = Ind(1:L_opt);         % The L_opt largest entries in S(:,k)
     
     %Ind_lamb = Ind_p(1:LAMBDA);   % THe first LAMBDA largest entries in S(:,k)
     
     s_hat = S(:,j)*0;  
     s_hat(Ind_p) = S(Ind_p,j); 
     
     e = Q*(S_opt(:,j)-s_hat); 
     snr = 20*log10(norm(Q*S_opt(:,j),2)/(norm(e,2)+1e-12)); 
     V_snr(j)=snr;
   % Checking the support 
     Ind_p = sort(Ind_p);          % descending order Ind_p as so is I_opt
     if norm(I_opt-Ind_p,1)==0     % Ind_p Different from I_opt
        Err_index(j)=1;            % successful !
     elseif snr > rho_snr          % srn > rho_snr dB
         Err_index(j)=1;           % successful ! 
     else
        Err_index(j)=0;            % failure !!!
     end

  %===== Check if the algprthm, yielding S, whether ===========
  %   ensures the indices of first LAMBDA lrgest entries fall in I_\s^\star

   %LAMBDA = 3;
%    E = intersect(Ind_lamb,I_opt);   % for checking if the first LAMBDA .. 
%    if length(E) == LAMBDA           % all in 
%       Err_LAMBDA(j)=1;              % the alg. works well 
%      else  
%       Err_LAMBDA(j)=0;              % NOT successful recovery 
%    end
                    
end       %   Output   s_0p, W_p
RecoverRate = sum(Err_index)/J;      % rate of successful recovery     
%RecoverRate_LAMBDA = sum(Err_LAMBDA)/J;
%================================================================
return 