function point = loong_y(x,y,z,x1,y1,z1,screen)
k = (screen - y1)/(y1-y);
x0 = (k+1) * x1 - k * x;
z0 = (k+1) * z1 - k * z;
point = [x0,z0];
end

%%
function point = loong_x(x,y,z,x1,y1,z1,screen)
k = (screen - x1)/(x1-x);
y0 = (k+1) * y1 - k * y;
z0 = (k+1) * z1 - k * z;
point = [y0,z0];
end

function point = loong_z(x,y,z,x1,y1,z1,screen)
k = (screen - z1)/(z1-z);
x0 = (k+1) * x1 - k * x;
y0 = (k+1) * y1 - k * y;
point = [x0,y0];
end