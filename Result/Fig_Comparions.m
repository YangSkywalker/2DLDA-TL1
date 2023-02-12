clear all; clc
figure(1)
% set parameter
x01 = -5:0.1:-0.02; x02 = 0.02:0.1:5;
x = -5:0.1:5;
y_TL1_100 = ((100+1).*abs(x))./(100+abs(x));
y_TL1_1 = ((1+1).*abs(x))./(1+abs(x));
y_TL1_01 = ((eps+1).*abs(x01))./(eps+abs(x01));
y_TL1_02 = ((eps+1).*abs(x02))./(eps+abs(x02));
y_TL1_0 = ((eps+1).*abs(x))./(eps+abs(x));

figure(1)
plot(x,y_TL1_0, 'color','g','LineStyle','-','linewidth',1.5); hold on
plot(x,y_TL1_1, 'color','r','LineStyle','-','linewidth',1.5); hold on
plot(x,y_TL1_100, 'color','b','LineStyle','-','linewidth',1.5); hold on
% plot(x02,y_TL1_02, 'color','g','LineStyle','-','linewidth',1.5); hold on
% plot(0, 0,'g.','markersize',12); hold on
set(gca,'XTicklabel',{'O'})
set(gca,'ytick','','xtick','');
set(gca,'xtick',0);
ylim([0,5])

set(gcf,'color','white'); %´°¿Ú±³¾°°×É«
backColor = [0.9843, 1.0000,0.9490];
set(gca, 'color', backColor); %´°¿Ú±³¾°°×É«
set(gca,'looseInset',[0 0 0 0]);
set(gca,'ticklength',[0 0]);
set(gca,'FontSize',12);
grid on;
ax = gca;
ax.GridColor = [0.7529    0.7529    0.7529];
ax.LineWidth = 1;
% text(-1.15,1.93,'asymptote y=a+1')
% legend({'L0','\mu_2(t)','\mu_1(t)','\mu_{p}(t)','\rho_a(t)'},'location','best')
legend({'a \rightarrow 0^{+}','a=1','a \rightarrow \infty'},'location', 'north','FontSize',8)

%%
clear all; clc
figure(2)
% set parameter
x = -5:0.1:5;
y2 = (abs(x)).^2;
y1 = (abs(x));
y_p = (abs(x)).^(0.5);
y_TL1 = ((1+1).*abs(x))./(1+abs(x));

plot(x,y_TL1, 'color','r','LineStyle','-','linewidth',1.5); hold on
plot(x,y2, 'color','g','LineStyle','-','linewidth',1.5); hold on
plot(x,y1, 'color','b','LineStyle','-','linewidth',1.5); hold on
plot(x,y_p, 'color','c','LineStyle','-','linewidth',1.5); hold on
plot(x,ones(size(x))*2,'color','r','LineStyle','--','linewidth',1)

set(gca,'ytick','','xtick','');
set(gca,'xtick',0);
set(gca,'XTicklabel',{'O'})
ylim([0,5])

set(gcf,'color','white'); %´°¿Ú±³¾°°×É«
backColor = [0.9843, 1.0000,0.9490];
set(gca, 'color', backColor); %´°¿Ú±³¾°°×É«
set(gca,'looseInset',[0 0 0 0]);
set(gca,'ticklength',[0 0]);
set(gca,'FontSize',12);
grid on;
ax = gca;
ax.GridColor = [0.7529    0.7529    0.7529];
ax.LineWidth = 1;

text(-1.15,1.93,'asymptote y=a+1')
legend({'TL_1-norm', 'L_2-norm','L_1-norm','L_p-norm'},'location', 'north','FontSize',8)

