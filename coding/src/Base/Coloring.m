function Color_matrix = Coloring(Color_matrix, pixel_idx, colorRGB, radius, fill_way)
if strcmp(fill_way, 'Square')
% Square Fill: Fills a square with an edge length of (2*radius-1) centered on the pixel
    rows = (pixel_idx(1) - (radius-1)) : (pixel_idx(1) + (radius-1));
    cols = (pixel_idx(3) - (radius-1)) : (pixel_idx(3) + (radius-1));

rows = rows(rows >= 1 & rows <= size(Color_matrix, 1)); 
cols = cols(cols >= 1 & cols <= size(Color_matrix, 2));  

if ~isempty(rows) && ~isempty(cols)
        color_3d = reshape(colorRGB, [1,1,3]);  
        Color_matrix(rows,cols, :) = repmat(color_3d, [length(rows), length(cols), 1]); 
end
    
elseif strcmp(fill_way, 'Round')
    % Round Fill: Fills a perfect circle centered on the pixel
        
        % 1. Calculate the valid region to avoid out-of-bound errors
        min_row = max(1, pixel_idx(1) - radius);
        max_row = min(size(Color_matrix, 1), pixel_idx(1) + radius);
        min_col = max(1, pixel_idx(3) - radius);
        max_col = min(size(Color_matrix, 2), pixel_idx(3) + radius);
        
        % 2. Generate grid only for the valid region
        [col_grid, row_grid] = meshgrid(min_col:max_col, min_row:max_row);
        
        % 3. Calculate squared distance from center
        distance_sq = (row_grid - pixel_idx(1)).^2 + (col_grid - pixel_idx(3)).^2;
        
        % 4. Create binary mask for the circle (1=inside, 0=outside)
        mask = distance_sq <= radius^2;
        col_grid = col_grid(mask);
        row_grid = row_grid(mask);
        color_3d = reshape(colorRGB, [1,1,3]);  
        Color_matrix(row_grid,col_grid, :) = repmat(color_3d, [length(row_grid), length(col_grid), 1]); 
else
    error('Unknown filling method, only "Square" or "Round" is supported');
end
end