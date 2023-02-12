function [W] = F2DLDA(X, Y, dim)  
%%
% % % % Input:
%   Data.X: Data matrix. Each d1xd2 matrix of Data.X is a sample.
%   Data.Y: Data label vector.
%   FunPara.MaxIter: Maximum iteration number.
%   FunPara.epsilon: Parameter epsilon > 0 is a small positive value.
%   FunPara.delta：10^-6
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

% X = Data.trainX; % 得到样本矩阵. Get the sample matrix.
% Y = Data.trainY; % 得到样本标签. Get the sample label.

[d1,d2,N] = size(X); % 得到矩阵的行数、列数和训练样本数. Obtain the number of rows and columns of the matrix and the number of training samples.

label = unique(Y); % 提取标签类别. Extract the label category.
c = length(label); % 计算训练样本类别数. Calculate the number of training sample categories.
N_i = zeros(c,1); % 用于存放每一类样本数量的矩阵. A matrix used to store the number of samples in each category.
X_mean = mean(X,3); % 计算训练样本均值. Calculate the mean value of training samples.

V = eye(d1,dim); % 初始化V；得到d1xdim的单位矩阵. Initialize V; The identity matrix of D1xdim is obtained.
I = eye(d1,d2); % 初始化I；得到d1xdim的单位矩阵. Initialize I; The identity matrix of D1xdim is obtained.
Sw = zeros(d1,d2); % 初始一个Sw；得到d1xd1的全零矩阵. Initial Sw; Get d1xd1 all zero matrix.

%%
for i = 1:c
    X_i_Index = find(Y==i); % 第i类训练样本的索引. Index of class i training samples.
    N_i(i) = length(X_i_Index); % 第i类训练样本的个数. Number of class i training samples.
    X_i = X(:,:,X_i_Index); % 第i个类对应的X数据. X data corresponding to the ith class.
    X_i_mean(:,:,i) = mean(X_i,3); % 第i个类的均值. The mean of the ith class.
    for s = 1:N_i(i)
        X_i_s = X_i(:,:,s); % 第i类的第s个训练样本. The s th training sample of class i.
        S = (X_i_s - X_i_mean(:,:,i)) * (X_i_s - X_i_mean(:,:,i))';
        Sw = Sw + S;
    end
    clear X_i_Index X_i
end
Sw = Sw / N;
if det(Sw) == 0
    Sw = Sw + (epsilon * I); % 该矩阵不可逆，加上一个正则项. The matrix is irreversible, plus a regular term.
else
    Sw = Sw; % 该矩阵可逆. The matrix is invertible.
end

%%
Q = chol(Sw); % 得到Q矩阵. Get the Q matrix.
H = zeros(d1,dim); % 设置一个全零的d1*dim的矩阵. Set a matrix of all zero d1*dim.
I_v = eye(d1,dim); % 另外设置一个d1*dim的单位矩阵. Also set up a d1*dim identity matrix.
for i = 1:c
    M_i(:,:,i) = (Q^(-1))' * (X_i_mean(:,:,i) - X_mean); % 计算Mi，求Q的逆：inv(Q) 或者 Q^(-1). Compute Mi and find the inverse of Q: inv(Q) or Q^(-1).
end

%% 开始迭代
for t = 1:MaxIter
    for i = 1:c
        r_i = 1/(trace((V' * M_i(:,:,i))' * (V' * M_i(:,:,i))))^0.5; % 计算ri. Compute ri.
        H_i = N_i(i) * r_i * (M_i(:,:,i) * M_i(:,:,i)' * V); 
        H = H + H_i; % 计算H. Compute H.
    end
    [U, S, P] = svd(H); % 对矩阵H进行奇异值分解. Singular value decomposition is performed on the matrix H.
    V_best = U * I_v * P';
    for i = 1:c
        K_best(i) = N_i(i) * (trace((V_best' * M_i(:,:,i))' * (V_best' * M_i(:,:,i))))^0.5;
        K(i) = N_i(i) * (trace((V' * M_i(:,:,i))' * (V' * M_i(:,:,i))))^0.5;
    end
    Eps = abs(sum(K_best) - sum(K)) / sum(K_best); % 终止条件判断. Termination condition judgment.
    if Eps <= eps 
        V = V_best;
        break; 
    end
    V = V_best;
end 
W = Q^(-1) * V;

end





