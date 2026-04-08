% This is the implementation of the algorithms based o  Iteratively 
% Reweighted Least Squares (IRLS) by R. Chartrand and WoTao. Yin in the paper:
% Iteratively Reweighted Algorithms for Compressive Sensing, in Proc.  
% of the 33rd IEEE ICASSP, 2008, pp. 3869¨C3872.

%  Problem :       z_opt = Argmin_z ||s||_phi_3  
%              s.t.  s = s_0 + W*z  -  our proposed characterization
%   where   phi_3(t)=(t^2+epsil)^(q/2) is the appro.g func. used in Yin's
%                                                       CY2008

function [s,z] = IRLS_alt3(s_0,W,q)   %(s_opt,Q)    %(s_0,W)       

%   [N,L]=size(Q);  x = Q*s_opt;
%   [U,Sig,V]=svd(Q);  Sig=Sig(:,1:N);     % Q = U[Sig 0]V^T
%   T_0 = V(:,1:N)*inv(Sig)*U';   W = V(:,N+1:L);     
%   s_0 = T_0*x;

[~, L_c]= size(W);
% Initials
  %z_p = zeros(L_c,1);   
           %lsqminnorm(Q,x);  % s = min_y||x-Qy||^2_2  with min_y||y||_2        
  s = s_0;  s_p=s_0*0;
  J=9;     epsilon = 10;     
for p=1:J
     epsilon =  epsilon*0.1;    pp=0;
     while norm(s-s_p,2)/norm(s,2) > sqrt(epsilon)/100   
                                           %Yin's criterion 
     pp = p+1  ;                  % when iter. no > 100, break while
     w = (abs(s).^2 + epsilon).^(q/2-1);       
     
     s_p = s; %z_p=z; 
     
   % To solve min_z||D{1/2}_w*s_0 + D^{1/2}_w*W*z||^2_2
      z = -inv(W'*(w.*W))*W'*(w.*s_0);    % the updated iterate 
%        b = -W'*(w.*s_0);  A=W'*(w.*W);       % b = A*z solved using 
%        tol=1e-6;
%        z = conjgrad(A,b,tol);      %       conjugate gradient
       
       s = s_0 + W*z;                   % the updated s
       if pp > 4500
           break
       end
     end
end



 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% subplot(3,1,1), stem(s)
% xlabel('s_est')
% subplot(3,1,2), stem(s_opt)
% xlabel('s_{opt}')
%  subplot(3,1,3), stem(s_opt-s)
%  xlabel('The diff')

end

 