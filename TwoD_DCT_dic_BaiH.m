function DCT_dict = TwoD_DCT_dic_BaiH(bb, RR)
    % Generate an 2-dimensional DCT dictionary 
% seting the parameters
  %bb=8; % block size: N = bb x bb
  %RR=4; % redundancy factor used for setting dictionary's size
  L=RR*bb^2;  
      % expected number of atoms in the dictionary: N x L BUT the dim.
      % of the generated dic is N x Pn^2 with 
      %                           Pn = ceil(RR^1/2*bb) >= L. See below 

% DCT dictionary
  Pn=ceil(sqrt(L));
  DCT=zeros(bb,Pn);   % bb x Pn
  for k=0:1:Pn-1
    V=cos([0:1:bb-1]'*k*pi/Pn);
    if k>0, V=V-mean(V); end
    DCT(:,k+1)=V/norm(V);
  end
  DCT_dict=kron(DCT,DCT);  % of size N x L_f with N = bb x bb and Pn^2 >=L
  
end