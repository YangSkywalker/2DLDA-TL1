 function [W] = D2L2BLDA(X,Y)  
% %SW:����
%�����
c = length(unique(Y));
%d1���У�N������������
[d1,~,N] = size(X);
SW=zeros(d1,d1);
for i=1:c
    SubclassIndex = find(Y==i);
    type(i) = length(SubclassIndex);
    %��i�����Ӧ��X���ݣ�
    SubData=X(:,:,SubclassIndex);
    %��i����ľ�ֵ��
    meani(:,:,i)=mean(SubData,3);
    Swi=zeros(d1,d1);
     %��i���������㵽��i����ľ�ֵ�ľ���ͣ�
    for j=1:type(i)
        Swi=Swi+(SubData(:,:,j)-meani(:,:,i))'*(SubData(:,:,j)-meani(:,:,i));
    end
    SW=SW+Swi;
end
%%%%%%%%%SB��Delta�����
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
                                              
