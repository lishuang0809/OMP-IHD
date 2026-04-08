% This is to plot
                
function[J] = plot_N_varying(J)

load data_N_varying   % save Data_N_varying     R T_c V_N

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%save data_N_e_varyingXvv V_2 N_e_2

F1=figure;
set(F1,'PaperUnits', 'centimeters'); 
set(F1,'Position',[0 0 600 400]);


F11=subplot(2,1,1);
set(F11,'LineWidth',1,'FontName','Times New Roman','FontSize',12,'FontWeight','Bold');

hold on;
plot(V_N,R(:,1),'-v')                %  虚线
hold on;
plot(V_N,R(:,2),'-*')                %  虚线
hold on;
plot(V_N,R(:,3),'-x')            %  虚线
hold on;
plot(V_N,R(:,4),'-o')               %  虚线
hold on;
plot(V_N,R(:,5),'-+')                %  虚线

axis([384 896  0.5 1.05])  
xlabel('(a)')
% ylabel('x(t),  x_1(t)');

% 
F12=subplot(2,1,2);
set(F12,'LineWidth',1,'FontName','Times New Roman','FontSize',12,'FontWeight','Bold');

hold on;
plot(V_N,T_c(:,1),'-v')                %  虚线

hold on;
plot(V_N,T_c(:,2),'-*')                %  虚线hold on;

hold on;
plot(V_N,T_c(:,3),'-x')                %  虚线

hold on;
plot(V_N,T_c(:,4),'-o')                %  虚线hold on;

hold on;
plot(V_N,T_c(:,5),'-+')                %  虚线


axis([384 896  0 150])
 xlabel('(b)')
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