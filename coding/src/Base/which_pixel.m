function which_pixel = which_pixels_y(x,z,px,pz,pixel)
num_x = ceil((px - x(1,1))/pixel);
num_z = ceil((pz - z(1,1))/pixel);
which_pixel = [num_z,num_x];
end

%%
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