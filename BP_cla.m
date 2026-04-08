% To address
%           s_opt = Argmin_s ||s||_0   s.t.     x = Q s
%    via        min_{s} ||s||_1       s.t. x = Q s 
%    with classical BP method (LP): define
%    z_1(i) = s(i), if s(i)>0; z_1(i)=0,     if s(i)<= 0.
%    z_2(i) = 0,    if s(i)>0; z_2(i)=-s(i), if s(i)<= 0
%    Clearly,  s = z_1 - z_2  & x = [Q -Q] z s.t. z = [z_1; z_2] => 0
%   Method is based on [FR2013]

                
function [s] = BP_cla(x,Q)       

% Inputs £¨x,Q,s_opt) can be generated using   
%                                        Rice_data_Generating.m

[N,L]=size(Q);


% Preparation for using linprog.m  - linear programming 
  f = ones(2*L,1);  
  A = - eye(2*L);  b = zeros(2*L,1);  % for A*z <= b  
  A_eq = [Q -Q];   b_eq = x;          % for A_eq*z = b_eq

% Choice I - slow and more memory taken but more accurate !
 
%   options = optimoptions('linprog','Algorithm','interior-point');
%   x_l=-ones(L_c+L_p,1)*100; x_u=-x_l; 
%   z = linprog(f,A,b,A_e,b_e,x_l,x_u,options);
     
% Choice II - fast and less meomory taken but less accurate  
  z = linprog(f,A,b,A_eq,b_eq);  % Linear Programming  to solve
                                 %  min f^T w  s.t.  A w <= b 
                                 % Here, w =[x; y] - (L_c + L_p) x 1, 
                                 %    f^Tw = sum(y)
  
 z_1=z(1:L); z_2=z(L+1:2*L);     % The optimal             
                            
 s = z_1 - z_2;        % 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%  [~,Ind]=sort(-abs(s));   
%  Ind = Ind(1:kappa);      % The kappa largest entries in s
%  s_est = zeros(L,1);
%  s_est(Ind) = s(Ind);   % kappa-sparse estimate of s_opt
% 
% 
%  
%  % Checking the support estimated based on s in terms of kappa largest 
%  %      based on    x = Q(:,Ind)*s_f
%  Q_f = Q(:,Ind); 
%  s_p = pinv(Q_f)*x;   
%  s_f = zeros(L,1);
%  s_f(Ind) = s_p;    % kappa-sparse estimate of s_opt
%  
%  error = norm(s_est-s_opt,2);
 % error_f = norm(s_f-s_est,2),
                      % should be NIL if the detected support is correct.
 
 
% subplot(3,1,1), stem(s)
% xlabel('s_est')
% subplot(3,1,2), stem(s_opt)
% xlabel('s_{opt}')
%  subplot(3,1,3), stem(s_opt-s)
%  xlabel('The diff')

%subplot(2,1,2), plot(e)
%xlabel('actual spare represent. error')
return

 