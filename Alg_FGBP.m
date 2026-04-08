% This is a revised version of the original Alg_GBP.m 

function [s] = Alg_FGBP(W,V_1,s_0,kappa,N_iter)

  [L,Lc] = size(W);       Q_bar = V_1*V_1';      
 % ================================================================

% =============         The algorithm         ==================  

%=====   Initial s  ======================
  [~,s] = BP_alt(s_0,W);
  residual = s;   
%===========================================
[~,index] = sort(abs(s),'descend');    
K_p = fix(0.75*kappa);
index = index(1:K_p);

for k=1:N_iter
    %========= Cheching if s is already is a good one   =====
     [~,index_p] = sort(abs(s),'descend');    
     if max(abs(s(index_p(kappa+1)))) < 10^(-4)  
         break
     end 
         
  % Computing  residual 
    Omega = [W -Q_bar(:,index)];  
    
   [z_bar,residual] = BP_alt(s_0,Omega);   % using BP_alt.m to solve
                           % solve min_beta||s_0 + Omega*z_bar||_1
   s=s_0*0;  s(index) = z_bar(Lc+1:Lc+length(index));  
   s = residual + s;
        
   % Updating the index
    [~,index] = sort(abs(s),'descend');    
    index = index(1:fix(1.05*kappa));
end

%=================End of the algorithm  ========================

return;

