function [A,V_k,indx]=OMP_LG(D,X,kappa)
% ========================================================
% Sparse coding of a group of signals based on a given dictionary and specified number of
% atoms to use.     X=DA
% input arguments:  D - the dictionary
%                   X - the signals to represent
%                   kappa - the maximum number of coefficients to use in OMP
% output arguments: A - sparse coefficient matrix.
% ========================================================

[~,P]=size(X);  %(^SS^):there are P signals with n dimensions
[~,K]=size(D);  %(^SS^):the dictionary is n-by-K
L = kappa;

D_c = normc(D);     % column l_2-normalized    !!!!

V_k = zeros(1,P);  % for storing the no. of iterations taken for each sample
t_1 = clock;  %tic 
 
for k=1:1:P    %(^SS^):对信号依次进行稀疏编码
    a=[];       %(^SS^):初始化解为零向量
    x=X(:,k);   %(^SS^):get the k-th signal and give it to x
    residual=x; %(^SS^):initialize the residual
    indx=[];  %zeros(L,1);%(^SS^):L - the maximum number of coefficients to use in OMP，indx-the solution support即非零元的位置
    for j=1:1:L
        proj=D_c'*residual;  % D_c   NOT D   !!!
        [~,pos]=max(abs(proj));
        pos=pos(1);          %只需要这样一个位置就好
        indx = [indx pos];  %indx(j)=pos;
        a=pinv(D(:,indx(1:j)))*x; % SL: changed "D" to "D_c"  %pinv 求伪逆，a与D(:,indx(1:j))有相同的维数
        residual=x-D(:,indx(1:j))*a; % SL: changed "D" to "D_c"
        if norm(residual,2) < 1e-12
            break;
        end
    end
    V_k(k) = j;
    temp=zeros(K,1);    %(^SS^):字典D有K列，则解a有K行，temp为所求解A的列向量初始化后的零向量
    temp(indx(1:j))=a;  %(^SS^):a是temp的非零部分
    A(:,k)=sparse(temp);%(^SS^):用temp去覆盖A的第k列
end

t_2=clock;   %toc 
t = etime(t_2,t_1);    % in second

%  LAMBDA = 1; rho_snr = 40;
%  RecoverRate = ComputingRecoveryRate(A,S_opt,Q,rho_snr,LAMBDA)
return;

