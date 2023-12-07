function [] = figure_pinocchio_bar_FDSVA_SO(time_ABA_deriv_SO_avx_g,...
                    time_ABA_deriv_SO_AD_cas_cg_avx_g,...
                    time_ABA_deriv_SO_ana_cg_avx_g,...
                    time_ABA_deriv_SO_avx_c,...
                    time_ABA_deriv_SO_AD_cas_cg_avx_c,...
                    time_ABA_deriv_SO_ana_cg_avx_c)



%% ABA derivs

figure

set(gcf,'units','inches','position',[1,1,6,4])
line_style = '-.';
ylim([5,1e4])

LW = 2.5;
MS = 3;
FS = 15;
FS_legend = 10; % 14
FS_x = 15; % 18
FS_text = 14;% 12
y3 = [time_ABA_deriv_SO_avx_g;...
    time_ABA_deriv_SO_ana_cg_avx_g;...
    time_ABA_deriv_SO_AD_cas_cg_avx_g;...
    time_ABA_deriv_SO_avx_c;...
    time_ABA_deriv_SO_ana_cg_avx_c;...
    time_ABA_deriv_SO_AD_cas_cg_avx_c].';

b=bar(y3);
b(1).FaceColor = 'g';
b(2).FaceColor = 'r';
b(3).FaceColor = 'b';
b(4).FaceColor = [0.8 1 0.8];
b(5).FaceColor = [1 0.8 0.8];
b(6).FaceColor = [0.77 0.77 1];

set(gca,'FontSize',FS)
set(gca, 'XTickLabel', {'double pend (2)','UR_{3}(6)','HyQ (18)','ATLAS (36)','Talos (50)'},'FontSize',FS_x)

set(gca, 'YScale', 'log')

%% text on top of bars

xt = get(gca, 'XTick');
hT=[];              % placeholder for text object handles

% for i=1:numel(b)  % iterate over number of bar objects
%    hT=[hT text(b(i).XData+b(i).XOffset,b(i).YData,num2str(ceil(b(i).YData.'),'%.d'), ...
%                           'VerticalAlignment','bottom','horizontalalign','right','FontSize',10)];                     
% end

xlabel('DOF ($n$)','interpreter','latex')
ylabel('run-time $(\mu s)$','Interpreter', 'latex')
grid on
grid minor
ylim([0.01,1e5])
yticks([1e-2,1e0,1e2,1e4])
l1 = legend('FDSVA-SO $\;\;\;\;\;\;\;\,$: GCC','FDSVA-SO (cg) : GCC',...
   'FD SO AD (cg)$\;\,$: GCC',...
   'FDSVA-SO $\;\;\;\;\;\;\;\,$: Clang','FDSVA-SO (cg) : Clang',...
   'FD SO AD (cg)$\;\,$: Clang');
set(l1,'Interpreter', 'latex','Location','northwest','Fontsize',FS_legend, 'NumColumns',2)

end
