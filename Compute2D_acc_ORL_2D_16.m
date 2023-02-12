%% ===== ORL_32X32 =====
% ------------ Original --------------
clear all; clc

% load data
load('E:\MatLab2019a\work\2DLDA-TL1\data\ORL_32_2D.mat');
data_path = 'E:\MatLab2019a\work\2DLDA-TL1\data\ORL_32_2D.mat';

% set parameter
ratio = 6/10;           % the ratio of train and test
block_size = 16;

iterations = 10;        % repeat 10 times
pcomp = 30;             % the number of principal component
maxDim = 30;            % maximal dimension 
space = linspace(1,maxDim,pcomp); 
interval = 1;

% store classification accuracy of D2LDATL1(a=100)
D2LDATL1_acc100 = zeros(iterations,pcomp);
% store classification accuracy of D2LDATL1(a=50)
D2LDATL1_acc50 = zeros(iterations,pcomp);
% store classification accuracy of D2LDATL1(a=10)
D2LDATL1_acc10 = zeros(iterations,pcomp);
% store classification accuracy of D2LDATL1(a=1)
D2LDATL1_acc1 = zeros(iterations,pcomp);
% store classification accuracy of D2LDATL1(a=0.5)
D2LDATL1_acc05 = zeros(iterations,pcomp);
% store classification accuracy of D2LDATL1(a=0.1)
D2LDATL1_acc01 = zeros(iterations,pcomp);
% store classification accuracy of D2LDATL1(a=0.05)
D2LDATL1_acc005 = zeros(iterations,pcomp);
% store classification accuracy of D2LDATL1(a=0.01)
D2LDATL1_acc001 = zeros(iterations,pcomp);
% store classification accuracy of D2LDATL1(a=0.001)
D2LDATL1_acc0001 = zeros(iterations,pcomp);
% store classification accuracy of D2LDA
D2LDA_acc = zeros(iterations,pcomp);
% store classification accuracy of L12DLDA
L12DLDA_acc = zeros(iterations,pcomp);
% store classification accuracy of Lp2DLDA(p=0.1)
Lp2DLDA_acc01 = zeros(iterations,pcomp);
% store classification accuracy of Lp2DLDA(p=0.2)
Lp2DLDA_acc02 = zeros(iterations,pcomp);
% store classification accuracy of Lp2DLDA(p=0.3)
Lp2DLDA_acc03 = zeros(iterations,pcomp);
% store classification accuracy of Lp2DLDA(p=0.4)
Lp2DLDA_acc04 = zeros(iterations,pcomp);
% store classification accuracy of Lp2DLDA(p=0.5)
Lp2DLDA_acc05 = zeros(iterations,pcomp);
% store classification accuracy of Lp2DLDA(p=0.6)
Lp2DLDA_acc06 = zeros(iterations,pcomp);
% store classification accuracy of Lp2DLDA(p=0.7)
Lp2DLDA_acc07 = zeros(iterations,pcomp);
% store classification accuracy of Lp2DLDA(p=0.8)
Lp2DLDA_acc08 = zeros(iterations,pcomp);
% store classification accuracy of Lp2DLDA(p=0.9)
Lp2DLDA_acc09 = zeros(iterations,pcomp);
% store classification accuracy of Lp2DLDA(p=1)
Lp2DLDA_acc1 = zeros(iterations,pcomp);


% store classification accuracy of F-2DLDA
F2DLDA_acc = zeros(iterations,pcomp);
% store classification accuracy of 2DBLDA
D2BLDA_acc = zeros(iterations,pcomp);

fprintf('Processing...... \n');
for i = 1:iterations
    % prepare the dataset of train and test
    [trainIdx, testIdx] = randomSplit2D(data_path, ratio);   
    x_train = blocksaltpepperPollute2D(X(:,:,trainIdx),block_size);  x_test = X(:,:,testIdx);
    y_train = Y(trainIdx);                                 y_test = Y(testIdx);
  
    % compute principal components of 2DLDA 
    start_time = clock;
    FunPara.d = 32;
    W_lda = D2LDA(x_train, y_train, FunPara, maxDim);  
    for j = space
        acc = knn_classifier2D(W_lda(:, 1:j), x_train, y_train, x_test, y_test);
        D2LDA_acc(i, j/interval) = acc;
    end
    end_time = clock; 
    fprintf('This is the %d-th iteration of 2DLDA on ORL_32X32, the elapsed time is %f s \n',i,etime(end_time,start_time));

    % compute principal components of D2LDATL1(a=100) 
    start_time = clock;
    W_TL1 = D2LDATL1_orth(x_train, y_train, maxDim, 100);
    % acc
    for j = space
        acc = knn_classifier2D(W_TL1(:,1:j),x_train,y_train,x_test,y_test);
        D2LDATL1_acc100(i,j/interval) = acc;
    end
    end_time = clock;
    fprintf('This is the %d-th iteration of 2DLDA-TL1(a=100) on ORL_32X32, the elapsed time is %f s \n',i,etime(end_time,start_time));

    % compute principal components of D2LDATL1(a=50) 
    start_time = clock;
    W_TL1 = D2LDATL1_orth(x_train, y_train, maxDim,50);
    % acc
    for j = space
        acc = knn_classifier2D(W_TL1(:,1:j),x_train,y_train,x_test,y_test);
        D2LDATL1_acc50(i,j/interval) = acc;
    end
    end_time = clock;
    fprintf('This is the %d-th iteration of 2DLDA-TL1(a=50) on ORL_32X32, the elapsed time is %f s \n',i,etime(end_time,start_time));

    % compute principal components of D2LDATL1(a=10) 
    start_time = clock;
    W_TL1 = D2LDATL1_orth(x_train, y_train,maxDim,10);
    % acc
    for j = space
        acc = knn_classifier2D(W_TL1(:,1:j),x_train,y_train,x_test,y_test);
        D2LDATL1_acc10(i,j/interval) = acc;
    end
    end_time = clock;
    fprintf('This is the %d-th iteration of 2DLDA-TL1(a=10) on ORL_32X32, the elapsed time is %f s \n',i,etime(end_time,start_time));

    % compute principal components of D2LDATL1(a=1) 
    start_time = clock;
    W_TL1 = D2LDATL1_orth(x_train, y_train,maxDim,1);
    % acc
    for j = space
        acc = knn_classifier2D(W_TL1(:,1:j),x_train,y_train,x_test,y_test);
        D2LDATL1_acc1(i,j/interval) = acc;
    end
    end_time = clock;
    fprintf('This is the %d-th iteration of 2DLDA-TL1(a=1) on ORL_32X32, the elapsed time is %f s \n',i,etime(end_time,start_time));
    
    % compute principal components of D2LDATL1(a=0.5) 
    start_time = clock;
    W_TL1 = D2LDATL1_orth(x_train, y_train,maxDim,0.5);
    % acc
    for j = space
        acc = knn_classifier2D(W_TL1(:,1:j),x_train,y_train,x_test,y_test);
        D2LDATL1_acc05(i,j/interval) = acc;
    end
    end_time = clock;
    fprintf('This is the %d-th iteration of 2DLDA-TL1(a=0.5) on ORL_32X32, the elapsed time is %f s \n',i,etime(end_time,start_time));
  
    % compute principal components of D2LDATL1(a=0.1) 
    start_time = clock;
    W_TL1 = D2LDATL1_orth(x_train, y_train,maxDim,0.1);
    % acc
    for j = space
        acc = knn_classifier2D(W_TL1(:,1:j),x_train,y_train,x_test,y_test);
        D2LDATL1_acc01(i,j/interval) = acc;
    end
    end_time = clock;
    fprintf('This is the %d-th iteration of 2DLDA-TL1(a=0.1) on ORL_32X32, the elapsed time is %f s \n',i,etime(end_time,start_time));
    
    % compute principal components of D2LDATL1(a=0.05) 
    start_time = clock;
    W_TL1 = D2LDATL1_orth(x_train, y_train,maxDim,0.05);
    % acc
    for j = space
        acc = knn_classifier2D(W_TL1(:,1:j),x_train,y_train,x_test,y_test);
        D2LDATL1_acc005(i,j/interval) = acc;
    end
    end_time = clock;
    fprintf('This is the %d-th iteration of 2DLDA-TL1(a=0.05) on ORL_32X32, the elapsed time is %f s \n',i,etime(end_time,start_time));
    
    % compute principal components of D2LDATL1(a=0.01) 
    start_time = clock;
    W_TL1 = D2LDATL1_orth(x_train, y_train,maxDim,0.01);
    % acc
    for j = space
        acc = knn_classifier2D(W_TL1(:,1:j),x_train,y_train,x_test,y_test);
        D2LDATL1_acc001(i,j/interval) = acc;
    end
    end_time = clock;
    fprintf('This is the %d-th iteration of 2DLDA-TL1(a=0.01) on ORL_32X32, the elapsed time is %f s \n',i,etime(end_time,start_time));
    
    % compute principal components of D2LDATL1(a=0.001) 
    start_time = clock;
    W_TL1 = D2LDATL1_orth(x_train, y_train,maxDim,0.001);
    % acc
    for j = space
        acc = knn_classifier2D(W_TL1(:,1:j),x_train,y_train,x_test,y_test);
        D2LDATL1_acc0001(i,j/interval) = acc;
    end
    end_time = clock;
    fprintf('This is the %d-th iteration of 2DLDA-TL1(a=0.001) on ORL_32X32, the elapsed time is %f s \n',i,etime(end_time,start_time));
   
    % compute principal components of L12DLDA 
    start_time = clock;
    delta = 0.05;
    W_L1 = L12DLDA(x_train, y_train, maxDim, delta);
    % acc
    for j = space
        acc = knn_classifier2D(W_L1(:,1:j),x_train,y_train,x_test,y_test);
        L12DLDA_acc(i,j/interval) = acc;
    end
    end_time = clock;
    fprintf('This is the %d-th iteration of L12DLDA on ORL_32X32, the elapsed time is %f s \n',i,etime(end_time,start_time));
 
    
    
    
    % compute principal components of Lp2DLDA£¨p=0.1£© 
    start_time = clock;
    W_Lp = Lp2DLDA(x_train, y_train, maxDim, 0.1);
    % acc
    for j = space
        acc = knn_classifier2D(W_Lp(:,1:j),x_train,y_train,x_test,y_test);
        Lp2DLDA_acc01(i,j/interval) = acc;
    end
    end_time = clock;
    fprintf('This is the %d-th iteration of Lp2DLDA(p=0.1) on ORL_32X32, the elapsed time is %f s \n',i,etime(end_time,start_time));    
    
    % compute principal components of Lp2DLDA£¨p=0.2£© 
    start_time = clock;
    W_Lp = Lp2DLDA(x_train, y_train, maxDim, 0.2);
    % acc
    for j = space
        acc = knn_classifier2D(W_Lp(:,1:j),x_train,y_train,x_test,y_test);
        Lp2DLDA_acc02(i,j/interval) = acc;
    end
    end_time = clock;
    fprintf('This is the %d-th iteration of Lp2DLDA(p=0.2) on ORL_32X32, the elapsed time is %f s \n',i,etime(end_time,start_time));    
    
    % compute principal components of Lp2DLDA£¨p=0.3£© 
    start_time = clock;
    W_Lp = Lp2DLDA(x_train, y_train, maxDim, 0.3);
    % acc
    for j = space
        acc = knn_classifier2D(W_Lp(:,1:j),x_train,y_train,x_test,y_test);
        Lp2DLDA_acc03(i,j/interval) = acc;
    end
    end_time = clock;
    fprintf('This is the %d-th iteration of Lp2DLDA(p=0.3) on ORL_32X32, the elapsed time is %f s \n',i,etime(end_time,start_time));    
   
    % compute principal components of Lp2DLDA£¨p=0.4£© 
    start_time = clock;
    W_Lp = Lp2DLDA(x_train, y_train, maxDim, 0.4);
    % acc
    for j = space
        acc = knn_classifier2D(W_Lp(:,1:j),x_train,y_train,x_test,y_test);
        Lp2DLDA_acc04(i,j/interval) = acc;
    end
    end_time = clock;
    fprintf('This is the %d-th iteration of Lp2DLDA(p=0.4) on ORL_32X32, the elapsed time is %f s \n',i,etime(end_time,start_time));    

    % compute principal components of Lp2DLDA£¨p=0.5£© 
    start_time = clock;
    W_Lp = Lp2DLDA(x_train, y_train, maxDim, 0.5);
    % acc
    for j = space
        acc = knn_classifier2D(W_Lp(:,1:j),x_train,y_train,x_test,y_test);
        Lp2DLDA_acc05(i,j/interval) = acc;
    end
    end_time = clock;
    fprintf('This is the %d-th iteration of Lp2DLDA(p=0.5) on ORL_32X32, the elapsed time is %f s \n',i,etime(end_time,start_time));    
    
    % compute principal components of Lp2DLDA£¨p=0.6£© 
    start_time = clock;
    W_Lp = Lp2DLDA(x_train, y_train, maxDim, 0.6);
    % acc
    for j = space
        acc = knn_classifier2D(W_Lp(:,1:j),x_train,y_train,x_test,y_test);
        Lp2DLDA_acc06(i,j/interval) = acc;
    end
    end_time = clock;
    fprintf('This is the %d-th iteration of Lp2DLDA(p=0.6) on ORL_32X32, the elapsed time is %f s \n',i,etime(end_time,start_time));    
    
    % compute principal components of Lp2DLDA£¨p=0.7£© 
    start_time = clock;
    W_Lp = Lp2DLDA(x_train, y_train, maxDim, 0.7);
    % acc
    for j = space
        acc = knn_classifier2D(W_Lp(:,1:j),x_train,y_train,x_test,y_test);
        Lp2DLDA_acc07(i,j/interval) = acc;
    end
    end_time = clock;
    fprintf('This is the %d-th iteration of Lp2DLDA(p=0.7) on ORL_32X32, the elapsed time is %f s \n',i,etime(end_time,start_time));    
    
    % compute principal components of Lp2DLDA£¨p=0.8£© 
    start_time = clock;
    W_Lp = Lp2DLDA(x_train, y_train, maxDim, 0.8);
    % acc
    for j = space
        acc = knn_classifier2D(W_Lp(:,1:j),x_train,y_train,x_test,y_test);
        Lp2DLDA_acc08(i,j/interval) = acc;
    end
    end_time = clock;
    fprintf('This is the %d-th iteration of Lp2DLDA(p=0.8) on ORL_32X32, the elapsed time is %f s \n',i,etime(end_time,start_time));    
    
    % compute principal components of Lp2DLDA£¨p=0.9£© 
    start_time = clock;
    W_Lp = Lp2DLDA(x_train, y_train, maxDim, 0.9);
    % acc
    for j = space
        acc = knn_classifier2D(W_Lp(:,1:j),x_train,y_train,x_test,y_test);
        Lp2DLDA_acc09(i,j/interval) = acc;
    end
    end_time = clock;
    fprintf('This is the %d-th iteration of Lp2DLDA(p=0.9) on ORL_32X32, the elapsed time is %f s \n',i,etime(end_time,start_time));    
    
    % compute principal components of Lp2DLDA£¨p=1£© 
    start_time = clock;
    W_Lp = Lp2DLDA(x_train, y_train, maxDim, 1);
    % acc
    for j = space
        acc = knn_classifier2D(W_Lp(:,1:j),x_train,y_train,x_test,y_test);
        Lp2DLDA_acc1(i,j/interval) = acc;
    end
    end_time = clock;
    fprintf('This is the %d-th iteration of Lp2DLDA(p=1) on ORL_32X32, the elapsed time is %f s \n',i,etime(end_time,start_time));    

    % compute principal components of 2DBLDA 
    start_time = clock;
    W_lda = D2L2BLDA(x_train, y_train); W_lda = W_lda(:, 1:maxDim);
    for j = space
        acc = knn_classifier2D(W_lda(:, 1:j), x_train, y_train, x_test, y_test);
        D2BLDA_acc(i, j/interval) = acc;
    end
    end_time = clock; 
    fprintf('This is the %d-th iteration of 2DBLDA on ORL_32X32, the elapsed time is %f s \n',i,etime(end_time,start_time));
    
    % compute principal components of F-2DLDA 
    start_time = clock;
    W_lda = F2DLDA(x_train, y_train, maxDim);
    for j = space
        acc = knn_classifier2D(W_lda(:, 1:j), x_train, y_train, x_test, y_test);
        F2DLDA_acc(i, j/interval) = acc;
    end
    end_time = clock; 
    fprintf('This is the %d-th iteration of F2DLDA on ORL_32X32, the elapsed time is %f s \n',i,etime(end_time,start_time));
end 

save('C:\MatLab2019a\work\TL1LDA\Result\ORL_2D_12','D2LDATL1_acc100','D2LDATL1_acc50','D2LDATL1_acc10','D2LDATL1_acc1',...
    'D2LDATL1_acc05','D2LDATL1_acc01','D2LDATL1_acc005','D2LDATL1_acc001','D2LDATL1_acc0001','D2LDA_acc','L12DLDA_acc',...
    'Lp2DLDA_acc01','Lp2DLDA_acc02','Lp2DLDA_acc03','Lp2DLDA_acc04','Lp2DLDA_acc05','Lp2DLDA_acc06','Lp2DLDA_acc07',...
    'Lp2DLDA_acc08','Lp2DLDA_acc09','Lp2DLDA_acc1','F2DLDA_acc','D2BLDA_acc','space');
fprintf('Finshed \n');
