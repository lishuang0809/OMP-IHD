% This is the implementation of the algorithms based o  Iteratively 
% Reweighted Least Squares (IRLS) by R. Chartrand and WoTao. Yin in the paper:
% Iteratively Reweighted Algorithms for Compressive Sensing, in Proc.  
% of the 33rd IEEE ICASSP, 2008, pp. 3869¨C3872.

%  Problem:          s_opt = Argmin_s ||s||^q_q   s.t.  x = Q s
%             where    0 <= q < 2
                
function [s] = Yin_IRLS_cla3(x,Q,q)       
%               0  <= q < 2    should be held.
% Initials
  s = lsqminnorm(Q,x);        % s = min_y||x-Qy||^2_2  with min_y||y||_2
  epsilon = 1;         

  J=9;    % was 10
for p=1:J
    epsilon = epsilon*0.1; 
    if epsilon < 10^(-8)
        break
    end
    
   % The IRLS for an epsilon given above with s set before
    w = (abs(s).^2 + epsilon).^(q/2-1);  
    %w = (abs(s) + epsilon).^(q-2); 
    w = w.^(-1);
    s_p = s;                       % The previus iterate
    
    s = (w.*Q')*inv(Q*(w.*Q'))*x;    % the updated iterate 
    ct=0;
    while norm(s-s_p,2)/norm(s_p,2) > sqrt(epsilon)/100
        ct = ct+1;
      w = (abs(s).^2 + epsilon).^(q/2-1);  
       %w = (abs(s) + epsilon).^(q-2); 
       w = w.^(-1);
       s_p = s;
       s = (w.*Q')*inv(Q*(w.*Q'))*x; % the updated iterate 
       if ct>4500
           break
       end
    end
end



 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
%  
% subplot(3,1,1), stem(s)
% xlabel('s_est')
% subplot(3,1,2), stem(s_opt)
% xlabel('s_{opt}')
%  subplot(3,1,3), stem(s_opt-s)
%  xlabel('The diff')

return

 