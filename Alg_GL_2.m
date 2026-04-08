function [A,indx]=Alg_GL_2(Q_c,W_c,X,kappa)
% ========================================================
% This is to implement the original OMP_IHD - Eqn.(10) in 
% Section II-B ``an l_2-relaxation"

% Sparse coding of a group of signals based on a given dictionary and specified number of
% atoms to use.     X = Q_bar S       - See 
% input arguments:  Q_c = V_1',  W_c = \Sigma^{-1}U'
%                   X - the signals to represent
%                   kappa - the maximum number of coefficients to use in OMP
% output arguments: A - sparse coefficient matrix.
% ========================================================

[~,P]=size(X);  %(^SS^):there are P signals with n dimensions
[~,K]=size(Q_c);  %(^SS^):the dictionary is n-by-K
L = kappa;

X = W_c*X;            % 文中的 x_c       ！！！
                      % x_c = Q_c*s   see ``Implementation complexity"

for k=1:1:P           %(^SS^):对信号依次进行稀疏编码
    a=[];             %(^SS^):初始化解为零向量
    x=X(:,k);         %(^SS^):get the k-th signal and give it to x
    residual=x;       %(^SS^):initialize the residual
    indx=[];          % zeros(L,1);    %(^SS^):L - the maximum number of coefficients 
                             %to use in OMP，indx-the solution support即非零元的位置
    for j=1:1:L
        proj = Q_c'*residual;              % 
        [~,pos]=max(abs(proj));
        pos=pos(1);                        %   只需要这样一个位置就好
        indx = [indx pos];                 %   indx(j)=pos;
        
        a=pinv(Q_c(:,indx(1:j)))*x;
                                %pinv 求伪逆，a与D(:,indx(1:j))有相同的维数
        
        residual = x - Q_c(:,indx(1:j))*a;
        if norm(residual,2) < 1e-12
            break;
        end
    end
    temp=zeros(K,1);    %(^SS^):字典D有K列，则解a有K行，temp为所求解A的列向量初始化后的零向量
    temp(indx(1:j))=a;  %(^SS^):a是temp的非零部分
    A(:,k)=sparse(temp);%(^SS^):用temp去覆盖A的第k列
end
            % 

return;

