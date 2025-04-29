function screen = Exit_effect(screen_total,Color_matrix_origin,Occlusion_ratio_both_sides)

% params:
% screen_total: Screen grid array, here only using the x-coordinate grid.
% Color_matrix_origin: A color grid of the same size as the screen grid, corresponding to the screen pixel display color.
% Occlusion_ratio_both_sides: The value ranges from 0 to 1, representing the ratio of the black screen area to the display area for each screen.

x_total = screen_total{1};
Occlusion = round(Occlusion_ratio_both_sides * size(x_total,2));
Color_left_origin = zeros([size(x_total,1),Occlusion,3]);
Color_right_origin = zeros([size(x_total,1),Occlusion,3]);
Color_matrix = Color_matrix_origin(:,Occlusion + 1:end - Occlusion,:);
Color_matrix = cat(2,Color_left_origin,Color_matrix,Color_right_origin);
screen = Color_matrix;
end