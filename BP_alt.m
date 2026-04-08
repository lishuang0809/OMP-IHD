% To address
%                  z_opt = Argmin_z ||s_0 + W*z||_1   
                
function[z_opt,s] = BP_alt(s_0,W)          

                         % s = s_0 +W*z  
% Initializing

[L_p,L_c]=size(W);  

% To solve     min_x ||s_0 + W*z||_1 using  linprog.m
     [n1,n2]=size(s_0);  L_p = n1*n2;
     f=[zeros(L_c,1);ones(L_p,1)];  
     A=[-W -eye(L_p);W -eye(L_p)];   b=[s_0;-s_0];
    
   % Choice I - slow and more memory taken but more accurate !
%      options = optimoptions('linprog','Algorithm','interior-point');
%      A_e = zeros(1,L_c+L_p); b_e = 0; x_l=-ones(L_c+L_p,1)*100; x_u=-x_l; 
%      y = linprog(f,A,b,A_e,b_e,x_l,x_u,options);
     
   % Choice II - fast and less meomory taken but less accurate  
     y = linprog(f,A,b);  % Linear Programming  to solve
                          %  min f^T w  s.t.  A w <= b 
                          % Here, w =[z; y] - (L_c + L_p) x 1, 
                          %    f^Tw = sum(y)
     z_opt=y(1:L_c);      
     s = s_0 + W*z_opt;         
return

 