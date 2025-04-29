function Color_matrix_lines = Constructing_structure(eye,x_total,y_total,z_total,deep,way_1,way_2,way_3,way_4,way_5)

% params:
% eye: Viewpoint
% x_total,y_total,z_total: 3D screen mesh array
% deep : The depth of the skybox scene
% way_1 - way_5 : Skybox image path

addpath(genpath('../Base/'));
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
%×ó1ÓÒ2ÉÏ3ÏÂ4ÖÐ5
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
                point = Extension_line_intersection_x(eye(1),eye(2),eye(3),where(1),where(2),where(3),left_wall_x);
                pixel = which_pixels_x(y_pic1,z_pic1,point(1),point(2),1);
                try
                Color_matrix_lines(i,j,:) = left_wall_pic(pixel(1),pixel(2),:);
                catch ME
                end
            case 2
                point = Extension_line_intersection_x(eye(1),eye(2),eye(3),where(1),where(2),where(3),right_wall_x);
                pixel = which_pixels_x(y_pic2,z_pic2,point(1),point(2),1);
                try
                Color_matrix_lines(i,j,:) = right_wall_pic(pixel(1),pixel(2),:);
                catch ME
                end
            case 3
                point = Extension_line_intersection_z(eye(1),eye(2),eye(3),where(1),where(2),where(3),ceiling_z);
                pixel = which_pixels_z(x_pic3,y_pic3,point(1),point(2),1);
                try
                Color_matrix_lines(i,j,:) = ceiling_pic(pixel(1),pixel(2),:);
                catch ME
                end
            case 4
                point = Extension_line_intersection_z(eye(1),eye(2),eye(3),where(1),where(2),where(3),ground_z);
                pixel = which_pixels_z(x_pic4,y_pic4,point(1),point(2),1);
                try
                Color_matrix_lines(i,j,:) = ground_pic(pixel(1),pixel(2),:);
                catch ME
                end
            case 5
                point = Extension_line_intersection_y(eye(1),eye(2),eye(3),where(1),where(2),where(3),hiding_y);
                pixel = which_pixels_y(x_pic5,z_pic5,point(1),point(2),1);
                try
                Color_matrix_lines(i,j,:) = hiding_pic(pixel(1),pixel(2),:);
                catch ME
                end
        end
    end
end
end