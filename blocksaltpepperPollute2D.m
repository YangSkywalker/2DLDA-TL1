function pdata = blocksaltpepperPollute2D(data,num)
% function: add block noise to image data
% data -- image data; num -- control the size of block 

[m, n, N] = size(data);
pdata = [];

for i=1:N
    img = data(:,:,i);
    point_x = m - num + 1;
    point_y = n  - num + 1;
    noise_vector = sign(randn(num^2,1)); index = find(noise_vector<=0);
    noise_vector(index) = 0;
    % give the left-up point of block noise
    start_point(1) = randperm(point_x,1);
    start_point(2) = randperm(point_y,1);
    img(start_point(1):start_point(1)+num-1,start_point(2):start_point(2)+num-1)=reshape(noise_vector,num,num);
    pdata(:,:,i) = img;
end

