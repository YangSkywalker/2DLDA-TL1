function [W] = Lp2DLDA(X,Y,dim,p)  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lp2DLDA: Robust bilateral Lp-norm two-dimensional linear discriminant analysis
%
% useage: [W] = Lp2DLDA(X,Y,dim,p)
% 
% Input:
%    X: input of Data.
%    Y: the class label.
%    dim: the reduced dimension.
%    p: the selection of p in Lp-norm
% Output:
%    W: transformation matrix (left side).

% Reference:
%    Chun-Na Li, Yuan-Hai Shao, Zhen Wang, Nai-Yang Deng "Robust bilateral Lp-norm two-dimensional linear discriminant analysis" 
%    Submitted 2017
%
%    Version 1.2 --Aug/2017 
%
%    Written by Chun-Na Li (na1013na@163.com)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
itmax = 50;
[d,n,N] = size(X); % N samples£¬each sample is with d*n dimension.
c = length(unique(Y)); 
w = rand(d,1); % Random initialization.
wk = []; % The k-th projection vector.
W = []; % The final projection matrix.

for k = 1:dim  
    barX = mean(X,3);
    Xmean = zeros(d,n,c); 
    num = zeros(c+1,1);      
    Hi = zeros(d,n,c);
    Zij = zeros(d,n,c);
    for i = 1:c
        tempMatrix = X(:,:,Y==i);
        num(i+1,1) = size(tempMatrix,3);
        Hi(:,:,i) = sum(tempMatrix,3)/num(i+1,1)-barX;
        Xmean(:,:,i) = sum(tempMatrix,3)/num(i+1,1);
    end
    Zij = X - Xmean(:,:,Y(1:N));
    theta = rand*pi/2;
    obj0 = -1e-12;
    it = 0;
    while 1
        A = zeros(d,1);
        B = 0; 
        C = 0;
        D = zeros(d,1);
        G = 0;
        it = it +1;
        for i = 1:c
            Atemp =  cumsum(Hi(:,:,i)*(num(i+1,1)*diag(sign(w'*Hi(:,:,i)).*((abs(w'*Hi(:,:,i))).^(p-1)))),2);
            Atemp =  Atemp(:,n);
            A = A + Atemp;
            C = C + num(i+1,1)*(norm(w'*Hi(:,:,i),p)^p);
            for j = 1:num(i+1,1)
                B = B + norm(w'*Zij(:,:,num(i,1)+j),p)^p;
                Dtemp =  cumsum(Zij(:,:,num(i,1)+j)*diag(sign(w'*Zij(:,:,num(i,1)+j)).*((abs(w'*Zij(:,:,num(i,1)+j))).^(p-1))),2);
                Dtemp =  Dtemp(:,n);
                D = D + Dtemp;
                G = G + norm(w'*Zij(:,:,num(i,1)+j),p)^p;
            end
        end
        obj(it) = C/G;
        A = p*A;
        D = p*D;
        G = G^2;      
        grad = (A*B-C*D)/G;       
        gradproj = grad - (w'*grad)*w;
        gradproj = gradproj/norm(gradproj);
        wk = w*cos(theta) + gradproj*sin(theta); 
        if obj(it)>obj0
           theta = min(2*theta,pi/2);
        else
            theta = theta/2.0;
        end
        if norm(w-wk) < 1e-5 ||it>itmax
            break;
        end
        w = wk;
        obj0 = obj(it);        
    end
    for h = 1:N
        X(:,:,h) = X(:,:,h)-wk*wk'*X(:,:,h); 
    end
    W = [W,wk];
end
end
