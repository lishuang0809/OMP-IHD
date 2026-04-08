% This is to plot
                
function[kk] = plot_condi_Q(kk)

%load data_C_QA    % V_condi  R T
%load data_C_QC    % V_condi  R T
load data_C_QE    % V_condi  R T

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%save data_N_e_varyingXvv V_2 N_e_2
%subplot(2,1,1), 
%hold on;
if kk==0
semilogx(V_condi,R(:,1),'-v')                %  虚线
hold on;
semilogx(V_condi,R(:,2),'-*')                %  虚线
hold on;
semilogx(V_condi,R(:,3),'-x')            %  虚线
% hold on;
% semilogx(V_condi,R(:,4),'-o')               %  虚线
% hold on;
% semilogx(V_condi,R(:,5),'-+')                %  虚线

%axis([384 896  0.5 1.05])  
%xlabel('(a)')
%  saveas(gcf,"fig-N_varying")            % saved in fig-format
%  saveas(gcf,"fig-condi_Q_R",'epsc')     % saved in eps-format
else

% subplot(2,1,2), 
% hold on;
%t1 = log10(T(:,1)); t2 = log10(T(:,2)); t3 = log10(T(:,3)); 
%t4 = log10(T(:,4)); t5 = log10(T(:,5));

t1 = T(:,1); t2 = T(:,2); t3 = T(:,3); 

semilogx(V_condi,t1,'-v')                %  虚线
hold on;
semilogx(V_condi,t2,'-*')                %  虚线
hold on;
semilogx(V_condi,t3,'-x')            %  虚线
% hold on;
% semilogx(V_condi,t4,'-o')               %  虚线
% hold on;
% semilogx(V_condi,t5,'-+')                %  虚线

%axis([log10(5) 10^5  0.5 1.05])  
%xlabel('(b)')

%saveas(gcf,"fig-condi_Q_T",'epsc')     % saved in eps-format

end
% F1=figure;
% set(F1,'PaperUnits', 'centimeters'); 
% set(F1,'Position',[0 0 600 400]);
% 
% 
% F11=subplot(2,1,1);
% set(F11,'LineWidth',1,'FontName','Times New Roman','FontSize',12,'FontWeight','Bold');
% 
% hold on;
% semilogx(V_condi,R(:,1),'-v')                %  虚线
% hold on;
% semilogx(V_condi,R(:,2),'-*')                %  虚线
% hold on;
% semilogx(V_condi,R(:,3),'-x')            %  虚线
% hold on;
% semilogx(V_condi,R(:,4),'-o')               %  虚线
% hold on;
% semilogx(V_condi,R(:,5),'-+')                %  虚线
% 
% %axis([384 896  0.5 1.05])  
% xlabel('(a)')
% % ylabel('x(t),  x_1(t)');

% 

% ylabel('x(t),  x_1(t)');

% plot(V_L,V_R(1,:)/100, '-v')   %  铏氱嚎
% hold on;
% plot(V_L,V_R(2,:)/100, '-*')   % 鐐圭嚎
% hold on;
% plot(V_L,V_R(3,:)/100, '-x')   % 鐐瑰垝绾?
% hold on;
% plot(V_L,V_R(4,:)/100, '-o')   % 鐐圭嚎
% hold on;
% plot(V_L,V_R(5,:)/100, '-+')   % 鐐瑰垝绾?

%print('fig-K_varying','-deps')   %After the figure shown, it can be save as a

%  saveas(gcf,"fig-N_varying")            % saved in fig-format
%  saveas(gcf,"fig-N_varying",'epsc')     % saved in eps-format


end
%%%%%%