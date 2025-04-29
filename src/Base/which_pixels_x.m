function which_pixel = which_pixels_x(y,z,py,pz,pixel)
num_y = ceil((py - y(1,1))/pixel);
num_z = ceil((pz - z(1,1))/pixel);
which_pixel = [num_z,num_y];
end
