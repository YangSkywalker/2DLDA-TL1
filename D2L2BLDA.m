 function [W] = D2L2BLDA(X,Y)  
% %SW:类内
%类别数
c = length(unique(Y));
%d1：行；N：样本个数；
[d1,~,N] = size(X);
SW=zeros(d1,d1);
for i=1:c
    SubclassIndex = find(Y==i);
    type(i) = length(SubclassIndex);
    %第i个类对应的X数据：
    SubData=X(:,:,SubclassIndex);
    %第i个类的均值：
    meani(:,:,i)=mean(SubData,3);
    Swi=zeros(d1,d1);
     %第i个类样本点到第i个类的均值的距离和：
    for j=1:type(i)
        Swi=Swi+(SubData(:,:,j)-meani(:,:,i))'*(SubData(:,:,j)-meani(:,:,i));
    end
    SW=SW+Swi;
end
%%%%%%%%%SB及Delta：类间
SB=zeros(d1,d1);
DeltaB = 0; 
for i = 1:c
    for j = 1:c
        if i<j 
          DeltaB = DeltaB + sqrt(type(i)*type(j))/N*(meani(:,:,i)-meani(:,:,j))*(meani(:,:,i)-meani(:,:,j))';
          SB=SB+sqrt(type(i)*type(j))/N*(meani(:,:,i)-meani(:,:,j))*(meani(:,:,i)-meani(:,:,j))';
        end
    end
end
SB= -SB;
Delta = DeltaB/4;
% Delta = DeltaB;
S = SB+ Delta*SW;
S = -S;
S = (S+S')/2;
[W,~] = eig(S);                                
                                              
