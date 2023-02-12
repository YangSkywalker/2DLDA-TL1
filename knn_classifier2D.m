function Acc = knn_classifier2D(W,x_train,y_train,x_test,y_test)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% knn_classifier: find the nearest neighbor of x_test in x_train
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

train_num = size(x_train,3);
test_num = size(x_test,3);

% projection
x_train_prj = [];
for i = 1:train_num
    x_train_prj(:,:,i) = x_train(:,:,i)*W;
end
x_test_prj = [];
for i = 1:test_num
    x_test_prj(:,:,i) = x_test(:,:,i)*W;
end

label = zeros(test_num,1);

for i = 1:test_num
    temp = zeros(1, train_num);
    for j = 1:train_num
        temp(j) = norm(x_train_prj(:,:,j) - x_test_prj(:,:,i),'fro');
    end
    optindexes = find(min(temp) == temp); optindex = optindexes(1);
    label(i) = y_train(optindex);
end
Acc = sum(label == y_test)/numel(y_test);

