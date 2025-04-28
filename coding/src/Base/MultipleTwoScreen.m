function [screens,colors,Key_points,x_y_z_into_three] = MultipleTwoScreen(l1,l2,h,theta0,theta1,pixel_wide,pixel_high)
see_wide_12 = pixel_wide * cosd(abs(theta0 - 90));
see_wide_23 = pixel_wide * cosd(abs(theta1 - 90));
first_axis_x = l1 * sind(theta0);
first_axis_y = l1 * cosd(theta0);
second_axis_x = first_axis_x + l2 * sind(theta1);
second_axis_y = first_axis_y + l2 * cosd(theta1);
[x1, z1] = meshgrid(0:see_wide_12:first_axis_x, 0:pixel_high:h);
[x2, z2] = meshgrid(first_axis_x:see_wide_23:second_axis_x, 0:pixel_high:h);
y1 = cotd(theta0) * abs(x1);
y2 = cotd(theta1) * abs(x2 - first_axis_x) + first_axis_y;

x_total = [x1 x2];
y_total = [y1 y2];
z_total = [z1 z2];
screens = {x_total,y_total,z_total};
colors = uint8(zeros([size(x_total), 3]));

Corner1 = [0,0,0];
Corner2 = [second_axis_x,second_axis_y,0];
Corner3 = [0,0,h];
Corner4 = [second_axis_x,second_axis_y,h];
Edge1 = [first_axis_x,first_axis_y,0];
Edge2 = [first_axis_x,first_axis_y,h];
Key_points = {Corner1,Corner2,Corner3,Corner4,Edge1,Edge2};
x_y_z_into_three = {x1,x2,y1,y2,z1,z2};
end