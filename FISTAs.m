% To address  
%                      min_{s} 1/2 ||x - Qs||^2_2 + lamb*||s||_1   
% using Iterative Shrinkage/Thresholding Algorithm £¨£©ISTA)-based methods.
% See Paper ¡°AMP-inspired deep learning network for sparse linear inverse 
%         problems", IEEE Trans. SP, vol. 65, n0.16, Aug. 2017, p. 4293 - 

function[s_FISTA] = FISTAs(x,Q,kappa,lamb,N_iter)   

% Inputs - 

[~,L]=size(Q);   

%%%% The ISTA algorithm
% V_error_ISTA = zeros(N_iter,1);      %
% 
% beta = 0.75/(norm(Q,2))^2;     % step-size£º  0 < beta < 1/||Q||^2_2  !!!
% 
% s = zeros(L,1);              % zero initial 
% for t=1:N_iter  
%              v = x - Q*s;         % iteration residual sequence 
%     V_error_ISTA(t) = norm(v,2);  
%     
%   % updating s with soft shringkage 
%     r = s + beta*Q'*v;
%     v = abs(r)-lamb;     ind = find(v < 0);  
%     
%     s = v; s(ind) = zeros(length(ind),1); % soft thresholding updating
%     s = sign(r).*s;                       % updated s (for t+1)
% end
% s_ISTA = s;   %norm(s-s_opt,2)/norm(s_opt,2)
%%%%%%%%%%%%%%%%%%%%%%% End of ISTA    

%%%% The FISTA algorithm
V_error_FISTA = zeros(N_iter,1);      %

beta = 0.75/(norm(Q,2))^2;     % step-size£º  0 < beta < 1/||Q||^2_2  !!!

s = zeros(L,1);    s_p=s;          % zero initial s=s_0, s_p=s_{-1}
for t=0:N_iter-1  
             v = x - Q*s;         % iteration residual sequence v_t
    V_error_FISTA(t+1) = norm(v,2);  
    
  % updating s with soft shringkage 
    r = s + beta*Q'*v + ((t-2)/(t+1))*(s-s_p);   % r_{t+1}
    s_p = s; 
    
    v = abs(r)-lamb;     ind = find(v < 0);  
    
    s = v; s(ind) = zeros(length(ind),1); % soft thresholding updating
    s = sign(r).*s;                       % updated s: s_{t+1}
    
        [~,index] = sort(abs(s),'descend'); 
    if abs(s(index(kappa+1))) < 1*10^(-4)  %10^(-3)  
       break
    end 
end
s_FISTA = s;  %norm(s-s_opt,2)/norm(s_opt,2)
%%%%%%%%%%%%%%%%%%%%%%% End of FISTA    


% F1=figure;
% axes('LineWidth',1,'FontName','Times New Roman','FontSize',12,'FontWeight','Bold');
% set(F1,'PaperUnits', 'centimeters'); 
% set(F1,'Position',[0 0 600 300]);
% hold on;
% %stem(t,s,'k','linewidth',1);
% %ylabel('Normalized Amplitude');
% %axis([0 15 -3.5 3]);   %[1.05,1.8,-1,1]);
% 
%   subplot(4,1,1), plot(V_error_ISTA)
%   xlabel('cost function - residual')
%   subplot(4,1,2), plot(s_opt)
%   xlabel('True s_{opt}')
%   subplot(4,1,3), plot(s_ISTA)
%   xlabel('Estimated s_{est}')
%   subplot(4,1,4), plot(s_opt-s_ISTA)
%   xlabel('s_{opt} - s_{est}')
%   
%   
% F2=figure;
% axes('LineWidth',1,'FontName','Times New Roman','FontSize',12,'FontWeight','Bold');
% set(F2,'PaperUnits', 'centimeters'); 
% set(F2,'Position',[0 0 600 300]);
% hold on;
% %stem(t,s,'k','linewidth',1);
% %ylabel('Normalized Amplitude');
% %axis([0 15 -3.5 3]);   %[1.05,1.8,-1,1]);
% 
%   subplot(4,1,1), plot(V_error_FISTA)
%   xlabel('cost function - residual')
%   subplot(4,1,2), plot(s_opt)
%   xlabel('True s_{opt}')
%   subplot(4,1,3), plot(s_FISTA)
%   xlabel('Estimated s_{est}')
%   subplot(4,1,4), plot(s_opt-s_FISTA)
%   xlabel('s_{opt} - s_{est}')
  
return

 