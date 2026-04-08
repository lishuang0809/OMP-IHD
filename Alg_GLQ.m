% This is a revised version of the original Alg_GBP.m 

function [s] = Alg_GLQ(W,V_1,s_0,kappa)

  [L,Lc] = size(W);       Q_bar = V_1*V_1';      
 % ================================================================

% =============         The algorithm         ==================  

%=====   Initial s  ======================
  q=1/2;
  [s,~] = IRLS_alt3(s_0,W,q);
  residual = s;   
%===========================================
index = [];
for k=1:fix(1.0*kappa)
    %========= Cheching if s is already is a good one   =====
     [~,index_p] = sort(abs(s),'descend');    
     if max(abs(s(index_p(kappa+1)))) < 10^(-4)  
         break
     end 
    % detecting the next atom - updating index
     [~,pos] = max(abs(residual));
     index = [index pos(1)];  
 
  % Computing  residual 
    Omega = [W -Q_bar(:,index)];    kappa_p=length(index);
         
   %[z_bar,residual] = BP_alt(s_0,Omega);   % using BP_alt.m to solve
                           % solve min_beta||s_0 + Omega*z_bar||_1
  q=1/2;
  [residual,z_bar] = IRLS_alt3(s_0,Omega,q);  
  
   s=s_0*0;  s(index) = z_bar(Lc+1:Lc+length(index));  
   s = residual + s;
end

%=================End of the algorithm  ========================

return;

