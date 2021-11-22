
function [drawRectangleImage] = drawRectangleFrame(image,windowLocation,windowSize)
[row,col] = size(image); % 输入图像尺寸
x = windowLocation(1);%矩形框位置坐标，其格式为[x,y]
y = windowLocation(2);
height = windowSize(1);%矩形框尺寸，其格式为[height,width]，即[高度,宽度]
width = windowSize(2);
if((x<=row && y<=col)&&(height<=row && width<=col))
    disp('矩形框合法！');
    LabelLineColor = 255;          % 标记线颜色
    drawRectangleImage = image;
    topMost = x-height;                  % 矩形框上边缘
    botMost = x;        % 矩形框下边缘
    lefMost = y-width;                  % 矩形框左边缘
    rigMost = y;        % 矩形框右边缘
    
    drawRectangleImage(topMost:botMost,lefMost,:) = LabelLineColor; % 左边框
    drawRectangleImage(topMost:botMost,rigMost,:) = LabelLineColor; % 右边框
    drawRectangleImage(topMost,lefMost:rigMost,:) = LabelLineColor; % 上边框
    drawRectangleImage(botMost,lefMost:rigMost,:) = LabelLineColor; % 下边框
else
    disp('Rectangle illegal');
end
