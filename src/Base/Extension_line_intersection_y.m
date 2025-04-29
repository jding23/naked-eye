function point = Extension_line_intersection_y(x,y,z,x1,y1,z1,screen)
% 
k = (screen - y1)/(y1-y);
x0 = (k+1) * x1 - k * x;
z0 = (k+1) * z1 - k * z;
point = [x0,z0];
end