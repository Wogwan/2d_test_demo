clear                  
close all;            
im1 = imread('controller_3D_Barrier.jpg');
im2 = imread('controller_3D_Barrier_Zoomin.jpg');

figure(1)
imshow(im1)
figure(2)
h_im2 = imshow(im2)

infoIm1 = imfinfo('controller_3D_Barrier.jpg')
infoIm2 = imfinfo('controller_3D_Barrier.jpg')

% 小尺寸彩色照片在风景图当中的居中处理
[m1,n1,l1] = size(im1);
[m2,n2,l2] = size(im2);
t = zeros(m1,n1,l1); 
t = uint8(t);
t((m1/2-m2/2+1):(m1/2+m2/2),(n1/2-n2/2+1):(n1/2+n2/2),:) = im2 ;%做居中处理
C = imadd(1*t,im1);%乘以0.5对中间小照片做透明处理

C((m1/2-m2/2+1):(m1/2+m2/2),(n1/2-n2/2+1):(n1/2+n2/2),:) =... 
C((m1/2-m2/2+1):(m1/2+m2/2),(n1/2-n2/2+1):(n1/2+n2/2),:)-...
im1((m1/2-m2/2+1):(m1/2+m2/2),(n1/2-n2/2+1):(n1/2+n2/2),:);  %对pic_1乘以0.5做补偿处理

figure(3)
imshow(C);
