% The codes below is to deal with x = Qs + e, where e = Gamma*alpha with
%      Gamma^T*Gamma = I_N_e  and N_e <= N, based on our proposed
%      characterization: s = T_0 (x-e) + W*z
%                          = T_0*x + [W - T_0*Gamma]*[z^T alpha^T]^T

function[S,T_c] = NewAlgs_4_Noise(Q,S_opt,X,kappa,W_bar,S_0,V_1,alg,Gamma,W,T_0)
[N,L]= size(Q);
[~,J] = size(S_opt);
%==============================================
if alg==104     %
    fprintf('Run Yin_IRLS\n')
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
if alg==107     %  the original Alg_GBP - very slow one!
    fprintf('Run Alg_{GL1}\n')
    tic,
    for j=1:J
        s = Alg_GBP(W_bar,V_1,S_0(:,j),kappa);
        S(:,j)=s;
    end
    T_c = toc;   % the time taken in seconds since tic
end
%============================================


%==============================================
if alg==108     %  the original Alg_GBP - very slow one!
    fprintf('Run Alg_{FGL1}\n')
    N_iter=5;
    tic,
    for j=1:J
        % s = Alg_GBP_rev(W,V_1,S_0(:,j),kappa);
        s = Alg_FGBP(W_bar,V_1,S_0(:,j),kappa,N_iter);
        S(:,j)=s;
    end
    T_c = toc;   % the time taken in seconds since tic
end
%============================================

%=================================================
if alg==109     % psi(t) = [|t|^2+epsilon]^{q/2} - l_q-norm
    fprintf('Run Alg_{GLQ}\n')
    tic,
    for j=1:J
        s = Alg_GLQ(W_bar,V_1,S_0(:,j),kappa);
        S(:,j) = s;
    end
    T_c=toc;
end
%===================================================================

%=================================================
if alg==110     % psi(t) = |t|/[|t} + epsil]
    fprintf('Run Alg_{FGLQ}\n')
    N_iter=2;
    tic,
    for j=1:J
        %s=Alg_FGDM2(W,V_1,S_0(:,j),kappa,N_iter);
        s = Alg_FGLQ(W_bar,V_1,S_0(:,j),kappa,N_iter);
        S(:,j) = s;
    end
    T_c=toc;
end
%===================================================================

%=================================================
if alg==111     % psi(t) = |t|/[|t} + epsil]
    fprintf('Run Projection + Alg_{FGLQ}\n')
    N_iter=2;
    Xtilde = (eye(N)-Gamma*Gamma')*X;
    Qtilde = (eye(N)-Gamma*Gamma')*Q;
    % Preparation
    [U,Sig,V]=svd(Qtilde);  Sig=Sig(:,1:N);        % Q = U[Sig 0]V^T
    Ntilde = sum(diag(Sig) > 1e-3);   % Qtilde may not be full row rank
    N = Ntilde;
    Sig=Sig(:,1:N);
    V_1 = V(:,1:N);                % Q_bar = V_1*V_1';
    T_0 = V(:,1:N)/Sig*U';         %SL: was V(:,1:N)*inv(Sig)*U';
    W = V(:,N+1:L);
    S_0 = T_0*Xtilde;

    tic,
    for j=1:J
        s = Alg_FGLQ(W,V_1,S_0(:,j),kappa,N_iter);
        S(:,j) = s;
    end
    T_c=toc;
end
%===================================================================



%===============The code below does not work==================================
if alg==112     % psi(t) = |t|/[|t} + epsil]
    fprintf('Run Alternating Minimization + Alg_{FGLQ}\n')
    N_iter=2;

    S0 = S_opt+ 0.01*randn(size(S_opt));  % Initialization
    S = S0;
    ct = 0;
    S_p = S/2;
    tic,
    while norm(S-S_p,'fro')/norm(S_p,'fro') > 1e-3
        ct = ct + 1
        % fix S, update alpha
        alpha = pinv(Gamma)*(X-Q*S);
        % fix alpha, update S with
        Xtilde = X - Gamma*alpha;
        % Preparation
        S_0 = T_0*Xtilde;
        S_p = S;
        for j=1:J
            s = Alg_FGLQ(W,V_1,S_0(:,j),kappa,N_iter);
            S(:,j) = s;
        end
        if ct>4500
            fprintf('Reached maximal number of iteration!')
            break
        end



    end


    T_c=toc;
end
%=====================The code above does not work==============================================





%=================================================
if alg==116     % psi(t) = |t|/[|t} + epsil]
    fprintf('Run Robust PCA\n')
    lambda = 0.1;    % Q = I: 0.1
    tol = 1e-7;
    maxIter = 1e3;
    tic,
    %[~, S, ~] = exact_alm_rpca(X, lambda, tol, maxIter);
    [~, S, ~] = exact_alm_rpca_Gamma(Gamma, X, lambda, tol, maxIter);
    T_c=toc;
end
%===================================================================




end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% F1=figure;
% set(F1,'PaperUnits', 'centimeters');
% set(F1,'Position',[0 0 600 400]);
%
%
% F11=subplot(2,1,1);
% set(F11,'LineWidth',1,'FontName','Times New Roman','FontSize',12,'FontWeight','Bold');
% hold on;
%
% plot(kappa,V_4(:,1), '-v')   %  虚线
% hold on;
% plot(kappa,V_7(:,1), '-*')   % 点线
% hold on;
% plot(kappa,V_8(:,1), '-x')   % 点划�?
%
% axis([1 32  0 1])
%  xlabel('(a)')
% % ylabel('x(t),  x_1(t)');
%
%
% F12=subplot(2,1,2);
% set(F12,'LineWidth',1,'FontName','Times New Roman','FontSize',12,'FontWeight','Bold');
% hold on;
%
% plot(kappa,V_4(:,2), '-v')   %  虚线
% hold on;
% plot(kappa,V_7(:,2), '-*')   % 点线
% hold on;
% plot(kappa,V_8(:,2), '-x')   % 点划�?
%
% axis([1 32  0 20])
%  xlabel('(b)')
% % ylabel('x(t),  x_1(t)');
%
% %print('fig-K_varying','-deps')   %After the figure shown, it can be save as a
%
% saveas(gcf,"fig-K_varying")            % saved in fig-format
% saveas(gcf,"fig-K_varying",'epsc')     % saved in eps-format
%

%%%%%%