clear;clc;
close all;

l = 100; theta0 = 60;theta1 = 90;theta2 = 120;
k = 2;
eye = [150,-350,50];
pixel_wide = 0.1;
pixel_high = 0.1;
l1 = 100;
l2 = 100;
l3 = 100;
h = 100;
% 重新添加路径（使用绝对路径）
proj_root = fileparts(mfilename('fullpath')); % 获取当前脚本所在目录
addpath(genpath(fullfile(proj_root, 'src')));
addpath(genpath(fullfile(proj_root, 'src', 'Base')));
[screens,colors,Key_points,xyz123] = MultipleThreeScreen(l1,l2,l3,h,theta0,theta1,theta2,pixel_wide,pixel_high);

x_total = screens{1};
y_total = screens{2};
z_total = screens{3};
%%
ptCloud = pcread('F:\bunny\bunny\reconstruction\bun_zipper.ply');
gras = grass(0,50,20);
M = 800*ptCloud.Location;
M(:, [2, 3]) = M(:, [3, 2]);
M(:,2) = -M(:,2) + 200;
M(:,1) = M(:,1) + 180;
M(:,3) = M(:,3) - 30;
cycle = load("cycle.mat").cylinder_points;
cycle(:,1) = cycle(:,1) + 100;
cycle(:,2) = cycle(:,2) + 125;
cycle(:,3) = cycle(:,3) + 50;

B = M;
%%
[jiao1,jiao2,jiao3,jiao4,bian1,bian2,bian3,bian4] = deal(Key_points{:});
%%
screen = zeros(0, 3);
color_screen = zeros(0, 3);
points_time = [];
deep = 300;
way_1 = "C:\Users\华为\Desktop\naked-eye\stars.jpg";
way_2 = "C:\Users\华为\Desktop\naked-eye\stars.jpg";
way_3 = "C:\Users\华为\Desktop\naked-eye\stars.jpg";
way_4 = "C:\Users\华为\Desktop\naked-eye\stars.jpg";
way_5 = "C:\Users\华为\Desktop\naked-eye\stars.jpg";


eye = [150 ,-300,60];
Color_matrix_origin = Constructing_structure(eye,x_total,y_total,z_total,theta0,theta1,theta2,pixel_wide,deep,way_1,way_2,way_3,way_4,way_5);
Color_matrix = Exit_effect(screens,Color_matrix_origin, 0.2);

points3D = M;
colors = getColorMatrixByZ(points3D);
numFrames = 100;
for i = 1:numFrames
condition_all = false(0,1);
k = i/100.0;
%%
layer_cross_point_1 = Cross(eye,jiao1,k);
layer_cross_point_2 = Cross(eye,jiao2,k);
layer_cross_point_3 = Cross(eye,jiao3,k);
layer_cross_point_4 = Cross(eye,jiao4,k);
layer_cross_point_5 = Cross(eye,bian1,k);
layer_cross_point_6 = Cross(eye,bian2,k);
layer_cross_point_7 = Cross(eye,bian3,k);
layer_cross_point_8 = Cross(eye,bian4,k);
layer_cross_Key_points = {layer_cross_point_1,layer_cross_point_2,layer_cross_point_3,layer_cross_point_4,layer_cross_point_5,layer_cross_point_6,layer_cross_point_7,layer_cross_point_8};
%%
point_to_layer_condition = Condition_Three(points3D,layer_cross_Key_points,theta0,theta1,theta2);

[row_idx,~] = find(point_to_layer_condition == 1);

for m = 1:length(row_idx)

M0 = points3D(row_idx(m),:);

Pixel = (M0 + k*eye)/(1+k);
pixel_idx = which_pixel_Three(Pixel(1),Pixel(2),Pixel(3),pixel_wide,xyz123,theta0,theta1,theta2);

colorRGB = [255,105,180];
radius = 2;
Color_matrix = Coloring(Color_matrix, pixel_idx, colorRGB, radius, 'Round');
end

end


figure;
hold on;
scatter3(points3D(:,1),points3D(:,2),points3D(:,3), ...
    5,"red",'filled');
h1 = surf(x_total, y_total, z_total,'FaceAlpha', 0.8, 'CData',Color_matrix, 'FaceColor', 'texturemap');% 绘制表面，并获取句柄
shading interp;

hold off;

ax = gca;

% 设置摄像机位置
ax.CameraPosition = eye;  % 摄像机位于 [10, 10, 10]

% 设置摄像机目标点
ax.CameraTarget = [110,50,50];  % 摄像机看向原点 [0, 0, 0]

% 设置摄像机上方向
ax.CameraUpVector = [0, 0, 1];  % 上方向为 Z 轴

% 设置摄像机视角宽度
%%

figure;
hold on;
scatter3(points3D(:,1),points3D(:,2),points3D(:,3), ...
    5,points3D(:,3),'filled');
h1 = surf(x_total, y_total, z_total,'FaceAlpha', 0.8, 'CData',Color_matrix, 'FaceColor', 'texturemap');% 绘制表面，并获取句柄
shading interp;
hold off;
%%
function layer_cross_point = Cross(eye,pixel,k)
layer_cross_point = (k + 1)*pixel - k*eye;
end

function which_pixel = which_pixel_Three(px,py,pz,pixel_wide,xyz123,theta0,theta1,theta2)
see_wide_12 = pixel_wide * cosd(abs(theta0 - 90));
see_wide_23 = pixel_wide * cosd(abs(theta1 - 90));
see_wide_34 = pixel_wide * cosd(abs(theta2 - 90));
first_axis = xyz123{1}(1,end);
second_axis = xyz123{2}(1,end);
third_axis = xyz123{3}(1,end);
if px < first_axis
num_x = ceil((px - 0)/see_wide_12);
num_z = ceil((pz - 0)/pixel_wide);
elseif px < second_axis
num_x = ceil((px - first_axis)/see_wide_23) + ceil((first_axis)/see_wide_12);
num_z = ceil((pz - 0)/pixel_wide);
elseif px < third_axis
num_x = ceil((px - second_axis)/see_wide_34) + ceil((first_axis)/see_wide_12) + ceil((second_axis - first_axis)/see_wide_23);
num_z = ceil((pz - 0)/pixel_wide);
end
which_pixel = [num_z,num_x,num_x];
end

function downsampledPoints = randomDownsample(points, percentage)
    % points: N-by-3 matrix where each row represents the coordinates of a point
    % percentage: Downsampling ratio, e.g., 0.1 means retaining 10% of the points
    N = size(points, 1);
    numPointsToKeep = round(N * percentage);
    shuffledIndices = randperm(N);
    indicesToKeep = shuffledIndices(1:numPointsToKeep);
    downsampledPoints = points(indicesToKeep, :);
end

function [colorMatrix] = getColorMatrixByZ(points)
    % 输入参数：
    % points - 一个N x 3的矩阵，其中N是点的数量，每一行代表一个点的(x, y, z)坐标

    % 输出参数：
    % colorMatrix - 一个N x 3的矩阵，每一行代表对应点的颜色(RGB)

    % 获取点的数量
    numPoints = size(points, 1);

    % 初始化颜色矩阵
    colorMatrix = zeros(numPoints, 3);

    % 获取z轴坐标
    zCoordinates = points(:, 3);

    % 找到z轴坐标的最大值和最小值
    zMax = max(zCoordinates);
    zMin = min(zCoordinates);

    % 根据z轴坐标的大小来设置颜色矩阵
    for i = 1:numPoints
        % 归一化z轴坐标，映射到[0, 1]区间
        normalizedZ = (zCoordinates(i) - zMin) / (zMax - zMin);
        
        % 根据归一化的z轴坐标来设置颜色
        % z坐标较低时，红色分量较大；z坐标较高时，蓝色分量较大
        colorMatrix(i, 1) = 1 - normalizedZ; % 红色分量
        colorMatrix(i, 3) = normalizedZ;     % 蓝色分量
    end
end

function gras = grass(x,y,height)
numPoints = 10000; % 点的数量
xRange = 300;       % X轴范围
yRange = 300;       % Y轴范围
maxHeight = 10;   % 草的最大高度

% 生成随机点
x = xRange * rand(numPoints, 1); % X坐标
z = yRange * rand(numPoints, 1); % Y坐标
y = maxHeight * rand(numPoints, 1); % Z坐标（模拟草的高度）

% 创建点云对象
ptCloud = pointCloud([x, z, y]);

% 为点云设置绿色颜色
greenColor = uint8([0, 255, 0]); % RGB颜色值 [R, G, B]
colors = repmat(greenColor, numPoints, 1); % 为每个点分配绿色
ptCloud.Color = colors; % 将颜色添加到点云对象

gras = ptCloud;
end

function point = loong_y(x,y,z,x1,y1,z1,screen)
k = (screen - y1)/(y1-y);
x0 = (k+1) * x1 - k * x;
z0 = (k+1) * z1 - k * z;
point = [x0,z0];
end

function point = loong_x(x,y,z,x1,y1,z1,screen)
k = (screen - x1)/(x1-x);
y0 = (k+1) * y1 - k * y;
z0 = (k+1) * z1 - k * z;
point = [y0,z0];
end

function point = loong_z(x,y,z,x1,y1,z1,screen)
k = (screen - z1)/(z1-z);
    x0 = (k+1) * x1 - k * x;
y0 = (k+1) * y1 - k * y;
point = [x0,y0];
end

function which_pixel = which_pixels_y(x,z,px,pz,pixel)
num_x = ceil((px - x(1,1))/pixel);
num_z = ceil((pz - z(1,1))/pixel);
which_pixel = [num_z,num_x];
end

function which_pixel = which_pixels_x(y,z,py,pz,pixel)
num_y = ceil((py - y(1,1))/pixel);
num_z = ceil((pz - z(1,1))/pixel);
which_pixel = [num_z,num_y];
end

function which_pixel = which_pixels_z(x,y,px,py,pixel)
num_x = ceil((px - x(1,1))/pixel);
num_y = ceil((py - y(1,1))/pixel);
which_pixel = [num_y,num_x];
end

function Color_matrix_origin = background(way,time,eye,x1,x2,x3,y1,y2,y3,z1,z2,z3,theta0,theta1,theta2,pixel_wide)
videoReader = VideoReader(way);
image_frame = read(videoReader,5*time);
image_frame= flipud(image_frame);
image_frame = imresize(image_frame, 0.5);
[height,width,~] = size(image_frame);
see_wide_12 = pixel_wide * cosd(abs(theta0 - 90));
see_wide_23 = pixel_wide * cosd(abs(theta1 - 90));
see_wide_34 = pixel_wide * cosd(abs(theta2 - 90));
[x_pic,z_pic] = meshgrid(-width/2:1:width/2,-height/2:1:height/2);
y_pic = zeros(size(x_pic)) + 500;
color_first = uint8(zeros([size(x1), 3]));
color_second = uint8(zeros([size(x2), 3]));
color_third = uint8(zeros([size(x3), 3]));
[rows, cols] = size(x1);
for i = 1:rows
    for j = 1:cols
        x_corner = x1(i, j) + see_wide_12/2;
        z_corner = z1(i, j) + pixel_wide/2;
        y_corner = y1(i, j);
        point = loong_y(eye(1),eye(2),eye(3),x_corner,y_corner,z_corner,500);
        pixel = which_pixels_y(x_pic,z_pic,point(1),point(2),1);
        try
        color_first(i,j,:) = image_frame(pixel(1),pixel(2),:);
        catch
        end
    end
end

[rows, cols] = size(x2);
for i = 1:rows
    for j = 1:cols
        x_corner = x2(i, j) + see_wide_23/2;
        z_corner = z2(i, j) + pixel_wide/2; 
        y_corner = y2(i, j);
        point = loong_y(eye(1),eye(2),eye(3),x_corner,y_corner,z_corner,500);
        pixel = which_pixels_y(x_pic,z_pic,point(1),point(2),1);
        try
        color_second(i,j,:) = image_frame(pixel(1),pixel(2),:);
        catch 
        end
    end
end
[rows, cols] = size(x3);
for i = 1:rows
    for j = 1:cols
        x_corner = x3(i, j) + see_wide_34/2;
        z_corner = z3(i, j) + pixel_wide/2;
        y_corner = y3(i, j);
        point = loong_y(eye(1),eye(2),eye(3),x_corner,y_corner,z_corner,500);
        pixel = which_pixels_y(x_pic,z_pic,point(1),point(2),1);
        try
        color_third(i,j,:) = image_frame(pixel(1),pixel(2),:);
        catch
        end
    end
end
Color_matrix_origin = [color_first color_second color_third];
end


function Color_matrix_lines = Constructing_structure(eye,x_total,y_total,z_total,deep,way_1,way_2,way_3,way_4,way_5)
Color_matrix_lines = zeros(size(x_total,1),size(x_total,2),3,'uint8');
[rows, cols] = size(x_total);
mask = zeros(size(x_total));
left_wall_x = x_total(1,1);
right_wall_x = x_total(end,end);
ceiling_z = z_total(end,end);
ground_z = z_total(1,1);
hiding_y = deep;
left_wall_pic = imread(way_1);
left_wall_pic = flipud(left_wall_pic);
left_wall_pic = imresize(left_wall_pic, [ceiling_z - ground_z, deep]);
[height,width,~] = size(left_wall_pic);
[y_pic1,z_pic1] = meshgrid(0:1:width,0:1:height);
right_wall_pic = imread(way_2);
right_wall_pic = flipud(right_wall_pic);
right_wall_pic = imresize(right_wall_pic, [ceiling_z - ground_z, deep]);
[height,width,~] = size(right_wall_pic);
[y_pic2,z_pic2] = meshgrid(0:1:width,0:1:height);
ceiling_pic = imread(way_3);
ceiling_pic = flipud(ceiling_pic);
ceiling_pic = imresize(ceiling_pic, [deep,right_wall_x - left_wall_x]);
[height,width,~] = size(ceiling_pic);
[x_pic3,y_pic3] = meshgrid(0:1:width,0:1:height);
ground_pic = imread(way_4);
ground_pic = flipud(ground_pic);
ground_pic = imresize(ground_pic, [deep,right_wall_x - left_wall_x]);
[height,width,~] = size(ground_pic);
[x_pic4,y_pic4] = meshgrid(0:1:width,0:1:height);
hiding_pic = imread(way_5);
hiding_pic = flipud(hiding_pic);
hiding_pic = imresize(hiding_pic, [ceiling_z - ground_z,right_wall_x - left_wall_x]);
[height,width,~] = size(hiding_pic);
[x_pic5,z_pic5] = meshgrid(0:1:width,0:1:height);
%左1右2上3下4中5
for i = 1:rows
    for j=1:cols
        vector = [x_total(i,j),y_total(i,j),z_total(i,j)] - eye;
        if vector(1) < 0 & vector(3) >= 0
            t1 = (left_wall_x - eye(1))/vector(1);
            t2 = (hiding_y - eye(2))/vector(2);
            t3 = (ceiling_z - eye(3))/vector(3);
            [~,index] = min([t1, t2, t3]);
            switch index
                case 1
                    mask(i,j) = 1;
                case 2
                    mask(i,j) = 5;
                case 3
                    mask(i,j) = 3;
            end
        elseif vector(1) > 0 & vector(3) >= 0
            t1 = (right_wall_x - eye(1))/vector(1);
            t2 = (hiding_y - eye(2))/vector(2);
            t3 = (ceiling_z - eye(3))/vector(3);
            [~,index] = min([t1, t2, t3]);
            switch index
                case 1
                    mask(i,j) = 2;
                case 2
                    mask(i,j) = 5;
                case 3
                    mask(i,j) = 3;
            end
        elseif vector(1) < 0 & vector(3) < 0
            t1 = (left_wall_x - eye(1))/vector(1);
            t2 = (hiding_y - eye(2))/vector(2);
            t3 = (ground_z - eye(3))/vector(3);
            [~,index] = min([t1, t2, t3]);
            switch index
                case 1
                    mask(i,j) = 1;
                case 2
                    mask(i,j) = 5;
                case 3
                    mask(i,j) = 4;
           end
        elseif vector(1) > 0 & vector(3) < 0
            t1 = (right_wall_x - eye(1))/vector(1);
            t2 = (hiding_y - eye(2))/vector(2);
            t3 = (ground_z - eye(3))/vector(3);
            [~,index] = min([t1, t2, t3]);
            switch index
                case 1
                    mask(i,j) = 2;
                case 2
                    mask(i,j) = 5;
                case 3
                    mask(i,j) = 4;
            end            
        end
    end
end
[rows, cols] = size(mask);
for i=1:rows
    for j=1:cols
        where = [x_total(i,j),y_total(i,j),z_total(i,j)];
        switch mask(i,j)
            case 1
                point = loong_x(eye(1),eye(2),eye(3),where(1),where(2),where(3),left_wall_x);
                pixel = which_pixels_x(y_pic1,z_pic1,point(1),point(2),1);
                try
                Color_matrix_lines(i,j,:) = left_wall_pic(pixel(1),pixel(2),:);
                catch ME
                end
            case 2
                point = loong_x(eye(1),eye(2),eye(3),where(1),where(2),where(3),right_wall_x);
                pixel = which_pixels_x(y_pic2,z_pic2,point(1),point(2),1);
                try
                Color_matrix_lines(i,j,:) = right_wall_pic(pixel(1),pixel(2),:);
                catch ME
                end
            case 3
                point = loong_z(eye(1),eye(2),eye(3),where(1),where(2),where(3),ceiling_z);
                pixel = which_pixels_z(x_pic3,y_pic3,point(1),point(2),1);
                try
                Color_matrix_lines(i,j,:) = ceiling_pic(pixel(1),pixel(2),:);
                catch ME
                end
            case 4
                point = loong_z(eye(1),eye(2),eye(3),where(1),where(2),where(3),ground_z);
                pixel = which_pixels_z(x_pic4,y_pic4,point(1),point(2),1);
                try
                Color_matrix_lines(i,j,:) = ground_pic(pixel(1),pixel(2),:);
                catch ME
                end
            case 5
                point = loong_y(eye(1),eye(2),eye(3),where(1),where(2),where(3),hiding_y);
                pixel = which_pixels_y(x_pic5,z_pic5,point(1),point(2),1);
                try
                Color_matrix_lines(i,j,:) = hiding_pic(pixel(1),pixel(2),:);
                catch ME
                end
        end
    end
end
end


function angle_point = rotatePointCloudAroundArbitraryAxis(P,axis,delta_angle)
rows = size(P,1);
rotation_matrix = [cos(delta_angle),sin(delta_angle);-sin(delta_angle),cos(delta_angle)];
for i=1:rows
    after =[P(i,1) - axis(1),P(i,2) - axis(2)] * rotation_matrix;
    angle_point(i,:) = [after(1),after(2),0] + [axis(1),axis(2),P(i,3)];
end
end

function [screens,colors,Key_points,x_y_z_into_three] = MultipleThreeScreen(l1,l2,l3,h,theta0,theta1,theta2,pixel_wide,pixel_high)
see_wide_12 = pixel_wide * cosd(abs(theta0 - 90));
see_wide_23 = pixel_wide * cosd(abs(theta1 - 90));
see_wide_34 = pixel_wide * cosd(abs(theta2 - 90));
first_axis_x = l1 * sind(theta0);
first_axis_y = l1 * cosd(theta0);
second_axis_x = first_axis_x + l2 * sind(theta1);
second_axis_y = first_axis_y + l2 * cosd(theta1);
third_axis_x = second_axis_x + l3 * sind(theta2);
third_axis_y = second_axis_y + l3 * cosd(theta2);
[x1, z1] = meshgrid(0:see_wide_12:first_axis_x, 0:pixel_high:h);
[x2, z2] = meshgrid(first_axis_x:see_wide_23:second_axis_x, 0:pixel_high:h);
[x3, z3] = meshgrid(second_axis_x:see_wide_34:third_axis_x, 0:pixel_high:h);
y1 = cotd(theta0) * abs(x1);
y2 = cotd(theta1) * abs(x2 - first_axis_x) + first_axis_y;
y3 = cotd(theta2) * abs(x3 - second_axis_x) + second_axis_y;

x_total = [x1 x2 x3];
y_total = [y1 y2 y3];
z_total = [z1 z2 z3];
screens = {x_total,y_total,z_total};
colors = uint8(zeros([size(x_total), 3]));

Corner1 = [0,0,0];
Corner2 = [third_axis_x,third_axis_y,0];
Corner3 = [0,0,h];
Corner4 = [third_axis_x,third_axis_y,h];
Edge1 = [first_axis_x,first_axis_y,0];
Edge2 = [second_axis_x,second_axis_y,0];
Edge3 = [first_axis_x,first_axis_y,h];
Edge4 = [second_axis_x,second_axis_y,h];
Key_points = {Corner1,Corner2,Corner3,Corner4,Edge1,Edge2,Edge3,Edge4};
x_y_z_into_three = {x1,x2,x3,y1,y2,y3,z1,z2,z3};
end

function screen = Exit_effect(screen_total,Color_matrix_origin,Occlusion_ratio_both_sides)
x_total = screen_total{1};
Occlusion = round(Occlusion_ratio_both_sides * size(x_total,2));
Color_left_origin = zeros([size(x_total,1),Occlusion,3]);
Color_right_origin = zeros([size(x_total,1),Occlusion,3]);
Color_matrix = Color_matrix_origin(:,Occlusion + 1:end - Occlusion,:);
Color_matrix = cat(2,Color_left_origin,Color_matrix,Color_right_origin);
screen = Color_matrix;
end

function Color_matrix = Coloring(Color_matrix, pixel_idx, colorRGB, radius, fill_way)
if strcmp(fill_way, 'Square')
% Square Fill: Fills a square with an edge length of (2*radius-1) centered on the pixel
    rows = (pixel_idx(1) - (radius-1)) : (pixel_idx(1) + (radius-1));
    cols = (pixel_idx(3) - (radius-1)) : (pixel_idx(3) + (radius-1));

rows = rows(rows >= 1 & rows <= size(Color_matrix, 1)); 
cols = cols(cols >= 1 & cols <= size(Color_matrix, 2));  

if ~isempty(rows) && ~isempty(cols)
        color_3d = reshape(colorRGB, [1,1,3]);  
        Color_matrix(rows,cols, :) = repmat(color_3d, [length(rows), length(cols), 1]); 
end
    
elseif strcmp(fill_way, 'Round')
    % Round Fill: Fills a perfect circle centered on the pixel
        
        % 1. Calculate the valid region to avoid out-of-bound errors
        min_row = max(1, pixel_idx(1) - radius);
        max_row = min(size(Color_matrix, 1), pixel_idx(1) + radius);
        min_col = max(1, pixel_idx(3) - radius);
        max_col = min(size(Color_matrix, 2), pixel_idx(3) + radius);
        
        % 2. Generate grid only for the valid region
        [col_grid, row_grid] = meshgrid(min_col:max_col, min_row:max_row);
        
        % 3. Calculate squared distance from center
        distance_sq = (row_grid - pixel_idx(1)).^2 + (col_grid - pixel_idx(3)).^2;
        
        % 4. Create binary mask for the circle (1=inside, 0=outside)
        mask = distance_sq <= radius^2;
        col_grid = col_grid(mask);
        row_grid = row_grid(mask);
        color_3d = reshape(colorRGB, [1,1,3]);  
        Color_matrix(row_grid,col_grid, :) = repmat(color_3d, [length(row_grid), length(col_grid), 1]); 
else
    error('Unknown filling method, only "Square" or "Round" is supported');
end
end

function point_to_layer_condition = Condition_Three(points3D,layer_cross_Key_points,theta0,theta1,theta2)
point_to_layer_condition = false(0,1);
point_change = points3D;
[layer_cross_point_1,layer_cross_point_2,layer_cross_point_3,layer_cross_point_4,layer_cross_point_5,layer_cross_point_6,layer_cross_point_7,layer_cross_point_8] = deal(layer_cross_Key_points{:});
for x = 1:size(point_change, 1)
    point_x = points3D(x, 1);
    point_y = points3D(x, 2);
    point_z = points3D(x, 3);

    if point_x < layer_cross_point_5(1) & point_x > layer_cross_point_1(1)  &  ((point_z > layer_cross_point_5(3))&(point_z < layer_cross_point_7(3)))
        point_change(x,1) = tand(90-theta0)*(point_x - layer_cross_point_1(1)) + layer_cross_point_1(2);
        condition = (point_change(x,1) < point_change(x,2)+2) & (point_change(x,1) > point_change(x,2)-2);
    elseif point_x < layer_cross_point_6(1)  &  ((point_z > layer_cross_point_5(3))&(point_z < layer_cross_point_7(3)))
        point_change(x,1) = tand(90-theta1)*(point_x - layer_cross_point_5(1)) + layer_cross_point_5(2);
        condition = (point_change(x,1) < point_change(x,2)+2) & (point_change(x,1) > point_change(x,2)-2);
    elseif point_x < layer_cross_point_2(1)  &  ((point_z > layer_cross_point_6(3))&(point_z < layer_cross_point_8(3)))
        point_change(x,1) = tand(90-theta2)*(point_x - layer_cross_point_2(1)) + layer_cross_point_2(2);
        condition = (point_change(x,1) < point_change(x,2)+2) & (point_change(x,1) > point_change(x,2)-2);
    else 
        condition = false;
    end
    point_to_layer_condition = [point_to_layer_condition;condition];
    
end
end