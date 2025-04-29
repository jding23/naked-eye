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