function f = TL12Dfun(x,y,w,a)
% function£ºcompute the objective function value 2DLDA-TL1 at w 
% Author: Skywalker Yang; Date: 2020/11/30

[m,n,N] = size(x);
label = unique(y);                 % class label
num_label = length(unique(y));     % the number of class label
f = 0;

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

% compute the numerator
numerator = 0;
for i = 1:num_label
    for k = 1:m
        numerator = numerator + num_eachclass(i) * ((a + 1) * abs((x_classmean(k,:,i) - x_mean(k, :)) * w) / ...
            (a + abs((x_classmean(k,:,i) - x_mean(k, :)) * w)));
    end
end

% compute the denominator
denominator = 0;
for i = 1:num_label
    for j = 1:length(class_label{1, i})
        for k = 1:m
            denominator = denominator + (a + 1) * abs((x(k,:,class_label{1, i}(j)) - x_classmean(k,:,i)) * w) / ....
                (a + abs((x(k,:,class_label{1, i}(j)) - x_classmean(k,:,i)) * w));
        end
    end
end

f = numerator / denominator;
end
