function [W] = D2LDA(X,Y,FunPara,dim)  
% =================================================
% % % % % % % % 2DLDA % % % % % %
% =================================================
c = length(unique(Y));
d = FunPara.d;
meanX=mean(X,3);%mean of matrices
SW=zeros(d,d);
SB=zeros(d,d);
for i=1:c
    SubclassIndex = find(Y==i);
    type(i) = length(SubclassIndex);
    SubData=X(:,:,SubclassIndex);
    meani(:,:,i)=mean(SubData,3);
    Swi=zeros(d,d);
    for j=1:type(i)
        Swi=Swi+(SubData(:,:,j)-meani(:,:,i))'*(SubData(:,:,j)-meani(:,:,i));
    end
    SW=SW+Swi;
    SB=SB+type(i)*(meani(:,:,i)-meanX)'*(meani(:,:,i)-meanX);
end
r=rank(SW);
% %     W 变换矩阵由v的最大的K-1个特征值所对应的特征向量构成
[W,D]=eig(SB,SW);
W=orth(W);
D = diag(D);
% % eigIdx1 = find(D < 1e-10);
% % eigIdx2 = find(D == Inf);
eigIdx1 = find(D < 1e-4);
eigIdx2 = find(D < 1e-4);
eigIdx = [eigIdx1;eigIdx2];
% eigIdx = eigIdx1;
D(eigIdx) = [];
W(:,eigIdx) = [];
[~, index] = sort(-D);
W = W(:,index);
if index>c
   index=1:c-1;
   W = W(:,index);
end
D = D(index);  
clear meanX SW SB SubclassIndex SubData Swi meani type;
%W = W(:,1:dim);?
 W = W(:,index);
 W = W(:,1:dim);
%=================================================================
% Requirement of energy is up to  FunPara.Ratio ( < 1),e.g. 0.90.
%=================================================================
%       FunPara.Ratio = 1;
%       sumEig = sum(D);
%       sumEig = sumEig*FunPara.Ratio;
%       sumNow = 0;
%       for idx = 1:length(D)
%           sumNow = sumNow + D(idx);
%           if sumNow >= sumEig
%               break;
%           end
%       end
%       W = W(:,1:idx); 
      
      
      
      