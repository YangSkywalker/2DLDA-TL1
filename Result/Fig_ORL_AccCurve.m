%% ORL Org
clear all; clc
% 2D
load('E:\MatLab2019a\work\2DLDA-TL1\Result\ORL_2D_Org.mat') 
space2D = space;
TL1_acc100_mean_2D = 100* mean(D2LDATL1_acc100); 
TL1_acc50_mean_2D = 100* mean(D2LDATL1_acc50);
TL1_acc10_mean_2D = 100* mean(D2LDATL1_acc10);
TL1_acc1_mean_2D = 100* mean(D2LDATL1_acc1);
TL1_acc05_mean_2D = 100* mean(D2LDATL1_acc05);
TL1_acc01_mean_2D = 100* mean(D2LDATL1_acc01);
TL1_acc005_mean_2D = 100* mean(D2LDATL1_acc005);
TL1_acc001_mean_2D = 100* mean(D2LDATL1_acc001);
TL1_acc0001_mean_2D = 100* mean(D2LDATL1_acc0001);

Lp2DLDA_acc01_2D = 100* mean(Lp2DLDA_acc01); 
Lp2DLDA_acc02_2D = 100* mean(Lp2DLDA_acc02);
Lp2DLDA_acc03_2D = 100* mean(Lp2DLDA_acc03);
Lp2DLDA_acc04_2D = 100* mean(Lp2DLDA_acc04);
Lp2DLDA_acc05_2D = 100* mean(Lp2DLDA_acc05);
Lp2DLDA_acc06_2D = 100* mean(Lp2DLDA_acc06);
Lp2DLDA_acc07_2D = 100* mean(Lp2DLDA_acc07);
Lp2DLDA_acc08_2D = 100* mean(Lp2DLDA_acc08);
Lp2DLDA_acc09_2D = 100* mean(Lp2DLDA_acc09);
Lp2DLDA_acc1_2D = 100* mean(Lp2DLDA_acc1);


D2LDA_acc_2D = 100* mean(D2LDA_acc);
L12DLDA_acc_2D = 100* mean(L12DLDA_acc);

F2DLDA_acc_2D = 100* mean(F2DLDA_acc);
D2BLDA_acc_2D = 100* mean(D2BLDA_acc);

figure(3)
plot(space2D,TL1_acc1_mean_2D,'r->','LineWidth',1.5,'MarkerSize',2.5); hold on
plot(space2D,Lp2DLDA_acc04_2D,'c-v','LineWidth',1.5,'MarkerSize',2.5); hold on
plot(space2D,D2LDA_acc_2D,'b-<','LineWidth',1.5,'MarkerSize',2.5); hold on
plot(space2D,L12DLDA_acc_2D,'m-^','LineWidth',1.5,'MarkerSize',2.5); hold on
plot(space2D,F2DLDA_acc_2D,'g-x','LineWidth',1.5,'MarkerSize',4); hold on
plot(space2D,D2BLDA_acc_2D,'-','color',[0.5412,0.1686,0.8863],'marker','d','LineWidth',1.5,'MarkerSize',2.5); hold on
xlabel('dimension'); ylabel('classification accuracy(%)');
ylim([35,100]);
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
legend({'2DLDA-TL1','Lp-2DLDA','2DLDA','L1-2DLDA','F2DLDA','2DBLDA'},'Location','southeast')

%% ORL 16X16
clear all; clc
% 2D
load('E:\MatLab2019a\work\2DLDA-TL1\Result\ORL_2D_16.mat') 
space2D = space;
TL1_acc100_mean_2D = 100* mean(D2LDATL1_acc100); 
TL1_acc50_mean_2D = 100* mean(D2LDATL1_acc50);
TL1_acc10_mean_2D = 100* mean(D2LDATL1_acc10);
TL1_acc1_mean_2D = 100* mean(D2LDATL1_acc1);
TL1_acc05_mean_2D = 100* mean(D2LDATL1_acc05);
TL1_acc01_mean_2D = 100* mean(D2LDATL1_acc01);
TL1_acc005_mean_2D = 100* mean(D2LDATL1_acc005);
TL1_acc001_mean_2D = 100* mean(D2LDATL1_acc001);
TL1_acc0001_mean_2D = 100* mean(D2LDATL1_acc0001);

Lp2DLDA_acc01_2D = 100* mean(Lp2DLDA_acc01); 
Lp2DLDA_acc02_2D = 100* mean(Lp2DLDA_acc02);
Lp2DLDA_acc03_2D = 100* mean(Lp2DLDA_acc03);
Lp2DLDA_acc04_2D = 100* mean(Lp2DLDA_acc04);
Lp2DLDA_acc05_2D = 100* mean(Lp2DLDA_acc05);
Lp2DLDA_acc06_2D = 100* mean(Lp2DLDA_acc06);
Lp2DLDA_acc07_2D = 100* mean(Lp2DLDA_acc07);
Lp2DLDA_acc08_2D = 100* mean(Lp2DLDA_acc08);
Lp2DLDA_acc09_2D = 100* mean(Lp2DLDA_acc09);
Lp2DLDA_acc1_2D = 100* mean(Lp2DLDA_acc1);

D2LDA_acc_2D = 100* mean(D2LDA_acc);
L12DLDA_acc_2D = 100* mean(L12DLDA_acc);

F2DLDA_acc_2D = 100* mean(F2DLDA_acc);
D2BLDA_acc_2D = 100* mean(D2BLDA_acc);

figure(5)
plot(space2D,TL1_acc1_mean_2D,'r->','LineWidth',1.5,'MarkerSize',2.5); hold on
plot(space2D,Lp2DLDA_acc09_2D,'c-v','LineWidth',1.5,'MarkerSize',2.5); hold on
plot(space2D,D2LDA_acc_2D,'b-<','LineWidth',1.5,'MarkerSize',2.5); hold on
plot(space2D,L12DLDA_acc_2D,'m-^','LineWidth',1.5,'MarkerSize',2.5); hold on
plot(space2D,F2DLDA_acc_2D,'g-x','LineWidth',1.5,'MarkerSize',4); hold on
plot(space2D,D2BLDA_acc_2D,'-','color',[0.5412,0.1686,0.8863],'marker','d','LineWidth',1.5,'MarkerSize',2.5); hold on
xlabel('dimension'); ylabel('classification accuracy(%)');
ylim([5,90]);
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
legend({'2DLDA-TL1','Lp-2DLDA','2DLDA','L1-2DLDA','F2DLDA','2DBLDA'},'Location','southeast')







