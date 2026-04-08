% To generate (v, Q, s_0) such that  v = Qs_0 with ||s_0||_0 <= kappa

% Conjugate gradient-related commands: 
%     See  bicg, bicgstab, bicgstabl, gmres, lsqr,  pcg, qmr, cgs 的文档

function [Q,S_opt,Q_condi]=gen_data_condQ(N,L,J,kappa,condi)
%=============================================  
 Q = randn(N,L);   Q = normc(Q);
              
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
  [U,S,V]=svd(Q);   
  rho = (1/condi)^(1/(N-1));    
  for n=1:N
      S(n,n)=rho^(n-1);
  end
  mu = (1-rho^(2*N))/(1-rho^2);
  c_1 = sqrt(L/mu);
  Q_condi = c_1*U*S*V';   % cond(Q_condi) = condi, ||Q_condi||^2_F=L
  %=========================================================
 
 
 
 return
 