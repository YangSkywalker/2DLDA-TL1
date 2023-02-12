function [W] = F2DLDA(X, Y, dim)  
%%
% % % % Input:
%   Data.X: Data matrix. Each d1xd2 matrix of Data.X is a sample.
%   Data.Y: Data label vector.
%   FunPara.MaxIter: Maximum iteration number.
%   FunPara.epsilon: Parameter epsilon > 0 is a small positive value.
%   FunPara.delta��10^-6
%   FunPara.eps: Tolerance eps : 10^-8.
%   dim: reduced dimension.
%
% % % % Eample:
%     Data.X = rand(32,32,10);
%     Data.Y = ones(10,1);
%     sFunPara.MaxIter = 50;
%     FunPara.epsilon = Autoregulation;
%     dim = 32;
%     [W] = F2DLDA(Data,FunPara,dim)
% 
% % % % Reference:
%    F2DLDA.
%    
%    Version 1.0 -- Nov/2020
%    Written by Yi-Fan Qi (yifanqicoco@163.com)

%%
eps = 10^-8; % Tolerance eps.
MaxIter = 50; % Maximum iteration number Itmax.
epsilon = 10^-6; % Where epsilon > 0 is a small positive value.

% eps = FunPara.eps; % Tolerance eps.
% MaxIter = FunPara.MaxIter; % Maximum iteration number Itmax.
% epsilon = FunPara.epsilon; % Where epsilon > 0 is a small positive value.

% X = Data.trainX; % �õ���������. Get the sample matrix.
% Y = Data.trainY; % �õ�������ǩ. Get the sample label.

[d1,d2,N] = size(X); % �õ������������������ѵ��������. Obtain the number of rows and columns of the matrix and the number of training samples.

label = unique(Y); % ��ȡ��ǩ���. Extract the label category.
c = length(label); % ����ѵ�����������. Calculate the number of training sample categories.
N_i = zeros(c,1); % ���ڴ��ÿһ�����������ľ���. A matrix used to store the number of samples in each category.
X_mean = mean(X,3); % ����ѵ��������ֵ. Calculate the mean value of training samples.

V = eye(d1,dim); % ��ʼ��V���õ�d1xdim�ĵ�λ����. Initialize V; The identity matrix of D1xdim is obtained.
I = eye(d1,d2); % ��ʼ��I���õ�d1xdim�ĵ�λ����. Initialize I; The identity matrix of D1xdim is obtained.
Sw = zeros(d1,d2); % ��ʼһ��Sw���õ�d1xd1��ȫ�����. Initial Sw; Get d1xd1 all zero matrix.

%%
for i = 1:c
    X_i_Index = find(Y==i); % ��i��ѵ������������. Index of class i training samples.
    N_i(i) = length(X_i_Index); % ��i��ѵ�������ĸ���. Number of class i training samples.
    X_i = X(:,:,X_i_Index); % ��i�����Ӧ��X����. X data corresponding to the ith class.
    X_i_mean(:,:,i) = mean(X_i,3); % ��i����ľ�ֵ. The mean of the ith class.
    for s = 1:N_i(i)
        X_i_s = X_i(:,:,s); % ��i��ĵ�s��ѵ������. The s th training sample of class i.
        S = (X_i_s - X_i_mean(:,:,i)) * (X_i_s - X_i_mean(:,:,i))';
        Sw = Sw + S;
    end
    clear X_i_Index X_i
end
Sw = Sw / N;
if det(Sw) == 0
    Sw = Sw + (epsilon * I); % �þ��󲻿��棬����һ��������. The matrix is irreversible, plus a regular term.
else
    Sw = Sw; % �þ������. The matrix is invertible.
end

%%
Q = chol(Sw); % �õ�Q����. Get the Q matrix.
H = zeros(d1,dim); % ����һ��ȫ���d1*dim�ľ���. Set a matrix of all zero d1*dim.
I_v = eye(d1,dim); % ��������һ��d1*dim�ĵ�λ����. Also set up a d1*dim identity matrix.
for i = 1:c
    M_i(:,:,i) = (Q^(-1))' * (X_i_mean(:,:,i) - X_mean); % ����Mi����Q���棺inv(Q) ���� Q^(-1). Compute Mi and find the inverse of Q: inv(Q) or Q^(-1).
end

%% ��ʼ����
for t = 1:MaxIter
    for i = 1:c
        r_i = 1/(trace((V' * M_i(:,:,i))' * (V' * M_i(:,:,i))))^0.5; % ����ri. Compute ri.
        H_i = N_i(i) * r_i * (M_i(:,:,i) * M_i(:,:,i)' * V); 
        H = H + H_i; % ����H. Compute H.
    end
    [U, S, P] = svd(H); % �Ծ���H��������ֵ�ֽ�. Singular value decomposition is performed on the matrix H.
    V_best = U * I_v * P';
    for i = 1:c
        K_best(i) = N_i(i) * (trace((V_best' * M_i(:,:,i))' * (V_best' * M_i(:,:,i))))^0.5;
        K(i) = N_i(i) * (trace((V' * M_i(:,:,i))' * (V' * M_i(:,:,i))))^0.5;
    end
    Eps = abs(sum(K_best) - sum(K)) / sum(K_best); % ��ֹ�����ж�. Termination condition judgment.
    if Eps <= eps 
        V = V_best;
        break; 
    end
    V = V_best;
end 
W = Q^(-1) * V;

end





