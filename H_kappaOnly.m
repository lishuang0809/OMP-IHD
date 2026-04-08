% This prog. is to  
%           i) find the kappa-sparse vector     s_hat = H_kappaOnly(s)
%                with the kappa largest  entriies in mangnitude  s

                
function[S_kappa,V_indx] = H_kappaOnly(S,kappa)          
            % S       - vectors to be checked
            % kappa   - sparssity level given  
            % S_kappa - This is the best estim. of kappa-sparse s*
% Initializing
                      [L,J]=size(S);  
  S_kappa = S*0;   
  V_indx = zeros(kappa,J);       % support of each S_kappa(:,j)
for k=1:J
   % The corresponding support of s to be examined
     [~,Ind]=sort(-abs(S(:,k)));   % sorting the entries 
     Ind_p = Ind(1:kappa);         % The kappa largest entries in S(:,k)
     S_kappa(Ind_p,k) = S(Ind_p,k);  
     V_indx(:,k) = sort(Ind_p);          % support of S_kappa
                             % This is the best estim. of kappa-sparse s*
end       %   Output   s_0p, W_p  
%================================================================
return 