 function W = L12DLDA(X,Y,dim,delta)  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L12DLDA: The L1-norm LDA for two-dimension redundency
%
% useage: [W] = L12DLDA(X,Y,itmax)
% 
% Input:
%    X: input of Data.
%    Y: the class label.
%    itmax: the iteration (No.) step.
% Output:
%    W: transfer matrix.
%
% Examples:
%    load('2Dexample.mat');
%    [W] = L12DLDA(X,Y)
% Reference:
%    Chun-Na Li, Yuan-Hai Shao, Nai-Yang Deng "Robust L1-norm two-dimensional 
%    linear discriminant analysis" Submitted 2014
%
%    Version 1.2 --Oct/2014 
%
%    Written by Chun-Na Li (na1013na@163.com)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (nargin <2 | nargin>10) % check correct number of arguments
    help L12DLDA
else
%     fprintf('_____________________________\n')
    if (nargin<3) itmax=1000; end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
itmax=100;
% X=X/255.0;
% XX = zeros(n,d,N);

%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % For right multiplication
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     for i = 1:size(X,3)
%         XX(:,:,i)  = X(:,:,i)';
%     end
%     X = XX;
%     clear XX;


[d,n,N]=size(X); % N samples£¬each sample is with d*n dimension.
c=size(unique(Y),1); % c classes.
w = rand(d,1); % Random initialization.
w = w/norm(w); % Normalize w.
wk=[]; % The k-th projection vector.
W=[]; % The final projection matrix.
% dim=30; % The maximum reduced dimension.
% delta=0.05; % The learning rate. 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following is to obtain W £¨for k=1:dim£©. The size of W is d times dim.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k=1:dim  
    %-------All mean and the mean of the i-th class------
    barX=zeros(d,n);
    barX=mean(X,3);% All mean.
    Xmean=zeros(d,n,c); % The matrix of means for c classes.
    num=zeros(c,1);      
     %-------N_i times Y_i------
    allNitimesYi=zeros(d,n,c);
    for i=1:c
        tempMatrix=X(:,:,find(Y==i)); % Put the samples of the i-th calss in tempMatrix.
        num(i,1)=size(tempMatrix,3);
        Xmean(:,:,i)=sum(tempMatrix,3)/num(i,1);% Compute the i-th mean and save it in Xmean(:,:,i).
        allNitimesYi(:,:,i)=num(i,1)*(Xmean(:,:,i)-barX);
    end
    %-------Z_ij------
    allZij=X-Xmean(:,:,Y(1:N));
    it=0;
    while 1
        it=it+1;
        b=zeros(d,1);
        p=zeros(d,1);
        r=ones(n,1);
        s=ones(n,1);
        numeratorw=0;
        denominatorw=0;
        %-------numerator of (6)------
        for i=1:c
            numeratorw=numeratorw+norm(w'*allNitimesYi(:,:,i),1);
        end
        %-------denominator of (6)------
        for h=1:N
            denominatorw=denominatorw+norm(w'*allZij(:,:,h),1);
        end
        %--------p(t)-------------
        for i=1:c
            temp=find(w'*(allNitimesYi(:,:,i))<0);
            r(temp')=-1;
            p=p+(allNitimesYi(:,:,i)*r);
        end
       %--------b(t)-------------
        for h=1:N
            temp=find(w'*(allZij(:,:,i))<0);
            s(temp')=-1;
            b=b+(allZij(:,:,h)*s);
        end
        %--------g(w(t))-------------
        if b==0
            g=p/(w'*p);
        else
            g=p/(w'*p)-b/(w'*b);
        end
        wk=w+delta*g;
        wk=wk/norm(wk);
        if wk'*b==0 || wk'*p==0
             wk=wk+(0.001+0.002*rand(d,1));
             wk=wk/norm(wk);
        end
        %--------compute new objective numerator and denominator of (6)-------------
        numeratorwk=0;
        denominatorwk=0;
        for i=1:c
           numeratorwk=numeratorwk+norm(wk'*allNitimesYi(:,:,i),1);
        end
        for h=1:N
           denominatorwk=denominatorwk+norm(wk'*allZij(:,:,h),1);
        end
        %-------convergence check-------
        if (abs(numeratorwk/denominatorwk-numeratorw/denominatorw)-(1e-6)<0)|| norm(w-wk) < 1e-6 ||it>itmax
            break;
        end
        w=wk;        
    end
     %-------Projcet samples in each recursive procedure------
    for h = 1:N
        X(:,:,h) = X(:,:,h)-wk*wk'*X(:,:,h);   % Makesure the projections are orthogonal to each other
    end
    W=[W,wk];
end