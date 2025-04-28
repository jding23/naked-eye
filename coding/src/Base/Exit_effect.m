function screen = Exit_effect(screen_total,Color_matrix,Occlusion_ratio_both_sides)
x_total = screen_total{1};
Occlusion = round(Occlusion_ratio_both_sides * size(x_total,2));
Color_left_origin = zeros([size(x_total,1),Occlusion,3]);
Color_right_origin = zeros([size(x_total,1),Occlusion,3]);
Color_matrix = Color_matrix([size(x_total,1),size(x_total,2) - 2*Occlusion ,3]);
Color_matrix = cat(2,Color_left_origin,Color_matrix,Color_right_origin);
screen = Color_matrix;
end