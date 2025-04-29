function which_pixel = which_pixels_z(x,y,px,py,pixel)
num_x = ceil((px - x(1,1))/pixel);
num_y = ceil((py - y(1,1))/pixel);
which_pixel = [num_y,num_x];
end