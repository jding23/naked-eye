function point_to_layer_condition = Condition_Two(points3D,layer_cross_Key_points,theta0,theta1)
point_to_layer_condition = false(0,1);
point_change = points3D;
[layer_cross_point_1,layer_cross_point_2,layer_cross_point_3,layer_cross_point_4,layer_cross_point_5,layer_cross_point_6] = deal(layer_cross_Key_points{:});
for x = 1:size(point_change, 1)
    point_x = points3D(x, 1);
    point_y = points3D(x, 2);
    point_z = points3D(x, 3);

    if point_x < layer_cross_point_5(1) & point_x > layer_cross_point_1(1)  &  ((point_z > layer_cross_point_5(3))&(point_z < layer_cross_point_6(3)))
        point_change(x,1) = tand(90-theta0)*(point_x - layer_cross_point_1(1)) + layer_cross_point_1(2);
        condition = (point_change(x,1) < point_change(x,2)+2) & (point_change(x,1) > point_change(x,2)-2);
    elseif point_x < layer_cross_point_4(1)  &  ((point_z > layer_cross_point_5(3))&(point_z < layer_cross_point_6(3)))
        point_change(x,1) = tand(90-theta1)*(point_x - layer_cross_point_5(1)) + layer_cross_point_5(2);
        condition = (point_change(x,1) < point_change(x,2)+2) & (point_change(x,1) > point_change(x,2)-2);
    else 
        condition = false;
    end
    point_to_layer_condition = [point_to_layer_condition;condition];
    
end
end