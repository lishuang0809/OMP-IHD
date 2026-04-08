% This prog. is to  evaluate the performance of algorithms    

function[S,T_c] = finalized_Master_4_Lis(X,Q,kappa,alg)    
      %             alg = {1, 2, 3, 4, 5, 6, 7}
      
       %[Q,S_opt]=generating_data(N,L,J,kappa);

% Preparation for (T_0, W) with Q for our charac.:  s = T_0*x + W*z
  [N,L]=size(Q);    [~,J] = size(X);
  
  [U,Sig,V]=svd(Q);  Sig=Sig(:,1:N);           % Q = U[Sig 0]V^T  
  V_1 = V(:,1:N);          % Q_bar = V_1*V_1';
  T_0 = V(:,1:N)/Sig*U';         %SL: was V(:,1:N)*inv(Sig)*U';   
  W = V(:,N+1:L);     
  W_c = Sig\U';                  %SL: was inv(Sig)*U';         
  %T_Q = T_0*Q;
  %=========================================================
                   S_0 = T_0*X;
                   
 %===========================================================

% Algorithms to be tested 
                         S=zeros(L,J); 
%====================================================================
if alg==1               % ==      classicl OMP  ==================
  tic,
    for j=1:J  
      s=OMP_LG(Q,X(:,j),kappa);
      S(:,j) = s;  
    end
  T_c = toc;   % the time taken in seconds since tic      
end
%==============================================================

%==============================================================
if alg==2    % ==      classic CoSaMP  ==================
  tic,
    for j=1:J  
        %s = BP_cla(X(:,j),Q)  ;  % SAME as  [~,s] = BP_alt(s_0,W); 
      s = CoSaMP(Q,X(:,j),kappa);
      S(:,j) = s;  
    end
  T_c = toc;   % the time taken in seconds since tic      
end
%====================================
  
%====================================
if alg==3                     % ==classicl BP_alt  ==================
  tic,
    for j=1:J  
      % [~,s] = BP_alt(S_0(:,j),W);  % SAME as  s = BP_cla(X(:,j),Q)  
       s = BP_cla(X(:,j),Q);     
      S(:,j) = s;  
    end
  T_c = toc;   % the time taken in seconds since tic      
 end
%====================================

%==============================================
if alg==4     % 
  tic,                       
     for j=1:J
                      q=1/2; 
        s = Yin_IRLS_cla3(X(:,j),Q,q)  ;
        S(:,j) = s;       % where  z_ini=0; kk=2; q=1/2;
     end
  T_c = toc;   % the time taken in seconds since tic      
end 
%===========================================
  
%==============================================
if alg==5     % 
  tic,                       
     for j=1:J
        lamb = 0.0001;   N_iter=2500;
        s = FISTAs(X(:,j),Q,kappa,lamb,N_iter);  
        S(:,j) = s;    
     end
  T_c = toc;   % the time taken in seconds since tic      
end 
%===========================================

%==================== Below are proposed algorithms

%==============================================
if alg==6     %  the proposed OMP_IHD applied to s_0 = Q_bar s
              % Named as Alg_{GL2} in the text
  tic,                       
       S = OMP_ILG(V_1*V_1',S_0,kappa);   % the origial OMP_IHD
       
                                  % Equiv. to OMP_ILG, but simplified
  T_c = toc;   % the time taken in seconds since tic      
end 
%===========================================  


%=========       Below are ours      =====================

%==============================================
if alg==7     %  the original Alg_GBP - very slow one! 
 tic,         % Named as Alg_{GL1} in the text
 for j=1:J
   s = Alg_GBP(W,V_1,S_0(:,j),kappa);
   S(:,j)=s;
 end
 T_c = toc;   % the time taken in seconds since tic   
end
%============================================

%==============================================
if alg==8     %  the original Alg_GBP - very slow one! 
   N_iter=5;   % Named as Alg_{GL1}^F in the text
 tic,   
 for j=1:J
   % s = Alg_GBP_rev(W,V_1,S_0(:,j),kappa);
   s = Alg_FGBP(W,V_1,S_0(:,j),kappa,N_iter);
   S(:,j)=s;
 end
 T_c = toc;   % the time taken in seconds since tic   
end
%============================================

%=================================================
if alg==9     % psi(t) = [|t|^2+epsilon]^{q/2} - l_q-norm

    tic,                       
    for j=1:J
       s = Alg_GLQ(W,V_1,S_0(:,j),kappa);
        S(:,j) = s;   
    end
    T_c=toc;   
end
%===================================================================
%=================================================
if alg==10     % psi(t) = |t|/[|t} + epsil]
          N_iter=2; 
    tic,                       
    for j=1:J
        %s=Alg_FGDM2(W,V_1,S_0(:,j),kappa,N_iter);
        s = Alg_FGLQ(W,V_1,S_0(:,j),kappa,N_iter);
        S(:,j) = s;   
    end
    T_c=toc;   
end
%===================================================================
end