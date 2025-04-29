function layer_cross_point = Cross_through_screens(eye,pixel,k)

% params:
% eye: Viewpoint.
% pixel:Screen grid (pixel position).
% k: k represents the ratio of the distance between the slice and the screen to the distance between the viewpoint and the screen, 
% which can determine the position of the slice.

layer_cross_point = (k + 1)*pixel - k*eye;

end