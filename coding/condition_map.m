function condition=condition_map(points_3D,layer_cross_points)
condition_all = false(0,1);
    layer_cross_point_1,
    layer_cross_point_2,
    layer_cross_point_3,
    layer_cross_point_4,
    layer_cross_point_5,
    layer_cross_point_6,
    layer_cross_point_7,
    layer_cross_point_8 = layer_cross_points;
    
for x = 1:size(points_3D, 1)
    point_x = points_3D(x, 1);
    point_y = points_3D(x, 2);
    point_z = points_3D(x, 3);

    if point_x < layer_cross_point_5(1) & point_x > layer_cross_point_1(1)  &  ((point_z > layer_cross_point_5(3))&(point_z < layer_cross_point_7(3)))
        points_3D(x,1) = tand(90-theta0)*(point_x - layer_cross_point_1(1)) + layer_cross_point_1(2);
        condition = (points_3D(x,1) < points_3D(x,2)+2) & (points_3D(x,1) > points_3D(x,2)-2);
    elseif point_x < layer_cross_point_6(1)  &  ((point_z > layer_cross_point_5(3))&(point_z < layer_cross_point_7(3)))
        points_3D(x,1) = tand(90-theta1)*(point_x - layer_cross_point_5(1)) + layer_cross_point_5(2);
        condition = (points_3D(x,1) < points_3D(x,2)+2) & (points_3D(x,1) > points_3D(x,2)-2);
    elseif point_x < layer_cross_point_2(1)  &  ((point_z > layer_cross_point_6(3))&(point_z < layer_cross_point_8(3)))
        points_3D(x,1) = tand(90-theta2)*(point_x - layer_cross_point_2(1)) + layer_cross_point_2(2);
        condition = (points_3D(x,1) < points_3D(x,2)+2) & (points_3D(x,1) > points_3D(x,2)-2);

    else 
        condition = false;
    end
    condition_all = [condition_all;condition];
    
end
condition = condition_all;
end