function point = Extension_line_intersection_z(x,y,z,x1,y1,z1,screen)
k = (screen - z1)/(z1-z);
x0 = (k+1) * x1 - k * x;
y0 = (k+1) * y1 - k * y;
point = [x0,y0];
end