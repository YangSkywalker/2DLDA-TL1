clear all; clc
figure(1)
load('E:\MatLab2019a\work\2DLDA-TL1\Result\Convergence2D.mat')
plot(Yale_2D_iter_org,Yale_2D_obj_org,'r-','linewidth',1.5); hold on
plot(Yale_2D_iter_noise,Yale_2D_obj_noise,'b-','linewidth',1.5); hold on

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
legend({'without noise' 'with noise'},'location', 'southeast');


%%
clear all; clc
figure(2)
load('E:\MatLab2019a\work\2DLDA-TL1\Result\Convergence2D.mat')
plot(ORL_2D_iter_org,ORL_2D_obj_org,'r-','linewidth',1.5); hold on
plot(ORL_2D_iter_noise,ORL_2D_obj_noise,'b-','linewidth',1.5); hold on
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
legend({'without noise' 'with noise'},'location', 'southeast');