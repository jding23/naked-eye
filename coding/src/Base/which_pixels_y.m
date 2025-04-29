function which_pixel = which_pixels_y(x,z,px,pz,pixel)
num_x = ceil((px - x(1,1))/pixel);
num_z = ceil((pz - z(1,1))/pixel);
which_pixel = [num_z,num_x];
end