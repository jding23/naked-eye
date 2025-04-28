function Color_matrix_origin = background(way,time,eye,x1,x2,x3,y1,y2,y3,z1,z2,z3,theta0,theta1,theta2,pixel_wide,deep)
videoReader = VideoReader(way);
image_frame = read(videoReader,5*time);
image_frame= flipud(image_frame);
image_frame = imresize(image_frame, 0.5);
[height,width,~] = size(image_frame);
see_wide_12 = pixel_wide * cosd(abs(theta0 - 90));
see_wide_23 = pixel_wide * cosd(abs(theta1 - 90));
see_wide_34 = pixel_wide * cosd(abs(theta2 - 90));
[x_pic,z_pic] = meshgrid(-width/2:1:width/2,-height/2:1:height/2);
y_pic = zeros(size(x_pic)) + deep;
color_first = uint8(zeros([size(x1), 3]));
color_second = uint8(zeros([size(x2), 3]));
color_third = uint8(zeros([size(x3), 3]));
[rows, cols] = size(x1);
for i = 1:rows
    for j = 1:cols
        x_corner = x1(i, j) + see_wide_12/2;
        z_corner = z1(i, j) + pixel_wide/2;
        y_corner = y1(i, j);
        point = loong_y(eye(1),eye(2),eye(3),x_corner,y_corner,z_corner,deep);
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
        point = loong_y(eye(1),eye(2),eye(3),x_corner,y_corner,z_corner,deep);
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
        point = loong_y(eye(1),eye(2),eye(3),x_corner,y_corner,z_corner,deep);
        pixel = which_pixels_y(x_pic,z_pic,point(1),point(2),1);
        try
        color_third(i,j,:) = image_frame(pixel(1),pixel(2),:);
        catch
        end
    end
end
Color_matrix_origin = [color_first color_second color_third];
end