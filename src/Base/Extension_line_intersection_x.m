function point = Extension_line_intersection_x(x,y,z,x1,y1,z1,screen)
k = (screen - x1)/(x1-x);
y0 = (k+1) * y1 - k * y;
z0 = (k+1) * z1 - k * z;
point = [y0,z0];
end
