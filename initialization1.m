function init_point=initialization1(x,y,a)

% function£ºfind the maximal obiective function value xi

[m,n]=size(x);
obj=zeros(n,1);

for i=1:n
    w = x(:,i)/norm(x(:,i));
    obj(i) = TL11Dfun(x,y,w,a);
end

[~,index]=max(obj);
init_point=x(:,index)/norm(x(:,index));
end
