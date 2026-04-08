% To generate (v, Q, s_0) such that  v = Qs_0 with ||s_0||_0 <= kappa

% Conjugate gradient-related commands: 
%     See  bicg, bicgstab, bicgstabl, gmres, lsqr,  pcg, qmr, cgs 的文档

function [Q,S_opt,T_0,W,mu_Q,kappa_h]=generating_data(N,L,J,kappa)
%=============================================  
 Q = randn(N,L);   G=Q'*Q; dd=diag(G); dd=sqrt(dd); dd=diag(dd);
 Q=Q/dd;    %   Q = Q*inv(dd)   - column-normalized
 %Q = eye(N);
 G=Q'*Q;  
              mu_Q=max(max(abs(G-eye(L)))); kappa_h=(1+1/mu_Q)/2;
              
 X = zeros(N,J);   % J samples of measurements v to be generated
 S_opt=zeros(L,J);   % J samples of sparse vectors s_opt
 
 for j=1:J
  s = zeros(L,1);
  loc = randperm(L,kappa);      %location vec. for the kappa non-zeros in s_0     
  s(loc) = randn(kappa,1);      %     % 
 % x=Q*s;     s = s/norm(x,2);   % s is normalized such that ||x||_2=1
 S_opt(:,j)=s;
 end
 
 % Preparation for (T_0, W) with Q for our charac.:  s = T_0*x + W*z
  [U,Sig,V]=svd(Q);  Sig=Sig(:,1:N);     % Q = U[Sig 0]V^T  
  V_1 = V(:,1:N);
  T_0 = V(:,1:N)*inv(Sig)*U';   W = V(:,N+1:L);     
  %T_Q = T_0*Q;
  %=========================================================
 
 
 
 return
 