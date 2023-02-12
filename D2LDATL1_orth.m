function [W] = D2LDATL1_orth(x, y, d, a)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% D2PCATL1: The TL1-norm LDA for two-dimensional data.
%
% Useage: [W] = D2LDATL1(x, y, d, a)
% 
% Input: x - matrix data
%        y - class label
%        d - number of projection vectors
%        a - a positive parameter(a>0)
%
% Output: W - transfer matrix
%
% Author: Skywalker Yang; Date: 2020/12/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[data_m,data_n,N] = size(x);  
label = unique(y);                 % class label
num_label = length(unique(y));     % the number of class label
W=[];
T = eye(data_n); 
dataX = x;
maxIter = 500;
tol = 1e-4;

for l = 1:d
    
    if(l == data_n)
         T = null(W'); w1=T;
         W=[W,w1];
         break;
    end
    
    % tranfer original data to low-dimensional subspace
    x = [];
    for i = 1:N
        x(:,:,i) = dataX(:,:,i) * T; 
    end
    [m,n,N] = size(x);
    % compute the mean value of every class and all samples
    class_label = cell(1, num_label);
    x_classmean = zeros(m,n,num_label);   % the mean value of every class
    x_mean = zeros(m,n);                  % the mean value of all samples
    num_eachclass = [];                   
    for i = 1:num_label
        class_label{1, i} = find(y == i); 
        num_eachclass(i) = length(class_label{1, i});
        for j = 1:length(class_label{1, i})
            x_classmean(:,:,i) = x_classmean(:,:,i) + x(:,:,class_label{1, i}(j));
            x_mean = x_mean + x(:,:,class_label{1, i}(j));
        end
        x_classmean(:,:,i) = x_classmean(:,:,i)/length(class_label{1, i});
    end
    x_mean = x_mean / N;
    
    w0 = rand(n,1); w0 = w0/norm(w0);
    theta0 = rand(1)*pi/2;
%     ('----%d-th------\n',l);
%     fprintf('%d\t%10.8f\n',0,TL12Dfun(x,y,w0,a));
    for kk = 1:maxIter
        %%% compute gradient
        % compute the numerator
        numerator1 = zeros(n,1);      % compute the the first term of the numerator 
        for i = 1:num_label
            for k = 1:m
                numerator1 = numerator1 + num_eachclass(i) * a * (a + 1) * sign((x_classmean(k,:,i) - x_mean(k, :)) * w0) * (x_classmean(k,:,i) - x_mean(k, :))' / ...
                    (a + abs((x_classmean(k,:,i) - x_mean(k, :)) * w0))^2;
            end
        end
        numerator2 = 0;      % compute the the second term of the numerator 
        for i = 1:num_label
            for j = 1:length(class_label{1, i})
                for k = 1:m
                    numerator2 = numerator2 + a * (a + 1) * abs((x(k,:,class_label{1, i}(j)) - x_classmean(k,:,i)) * w0) / ....
                        (a + abs((x(k,:,class_label{1, i}(j)) - x_classmean(k,:,i)) * w0));
                end
            end
        end
        numerator3 = zeros(n,1);      % compute the the third term of the numerator 
        for i = 1:num_label
            for j = 1:length(class_label{1, i})
                for k = 1:m
                    numerator3 = numerator3 + a * (a + 1) * sign((x(k,:,class_label{1, i}(j)) - x_classmean(k,:,i)) * w0) * (x(k,:,class_label{1, i}(j)) - x_classmean(k,:,i))' / ....
                        (a + abs((x(k,:,class_label{1, i}(j)) - x_classmean(k,:,i)) * w0))^2;
                end
            end
        end 
        numerator4 = 0;      % compute the the fourth term of the numerator 
        for i = 1:num_label
            for k = 1:m
                numerator4 = numerator4 + num_eachclass(i) * (a + 1) * abs((x_classmean(k,:,i) - x_mean(k, :)) * w0) / ...
                    (a + abs((x_classmean(k,:,i) - x_mean(k, :)) * w0));
            end
        end
        numerator = numerator1 * numerator2 - numerator3 * numerator4;
        % compute the denominator
        denominator = numerator2^2;
        
        % compute gradient
        gradf = numerator/denominator;
        
        g = gradf-(gradf'*w0)*w0;  g0=g/norm(g);
        while 1
            w1 = w0*cos(theta0) + g0*sin(theta0);
            if(TL12Dfun(x,y,w1,a)>=TL12Dfun(x,y,w0,a))
                break;
            end
            theta0=theta0/2;
        end
%         fprintf('%d\t%10.8f\n',kk,TL12Dfun(x,y,w1,a));
        %  update theta1
        theta1 = min(theta0*2,pi/2);
        if(TL12Dfun(x,y,w1,a) - TL12Dfun(x,y,w0,a) < tol)
            break;
        end
        theta0 = theta1;
        w0 = w1;
    end
    w1 = T * w1;
    W = [W, w1];
    T = null(W');
%     W(:,l) = w1;
%     for i = 1:N
%         x(:,:,i) = x(:,:,i) - x(:,:,i)*w1*w1';
%     end
end  
disp('2DLDA-TL1 finished');
end
