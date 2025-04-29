function which_pixel = which_pixel_Two(px,py,pz,pixel_wide,xyz123,theta0,theta1)
see_wide_12 = pixel_wide * cosd(abs(theta0 - 90));
see_wide_23 = pixel_wide * cosd(abs(theta1 - 90));
first_axis = xyz123{1}(1,end);
second_axis = xyz123{2}(1,end);
if px < first_axis
num_x = ceil((px - 0)/see_wide_12);
num_z = ceil((pz - 0)/pixel_wide);
elseif px < second_axis
num_x = ceil((px - first_axis)/see_wide_23) + ceil((first_axis)/see_wide_12);
num_z = ceil((pz - 0)/pixel_wide);
end
which_pixel = [num_z,num_x,num_x];
end