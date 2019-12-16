function Ih = apply_H(I, H) 
    
    % Apply H to corners of I to determine the horizontal and vertical 
    % limits of the output image
    [height, width, nChannels] = size(I);

    corner1 = H * [0; 0; 1];
    corner1_y = floor(corner1(1) / corner1(3));
    corner1_x = floor(corner1(2) / corner1(3));
    
    corner2 = H * [height; 0; 1];
    corner2_y = floor(corner2(1) / corner2(3));
    corner2_x = floor(corner2(2) / corner2(3));
    
    corner3 = H * [0; width; 1];
    corner3_y = floor(corner3(1) / corner3(3));
    corner3_x = floor(corner3(2) / corner3(3));
    
    corner4 = H * [height; width; 1];
    corner4_y = floor(corner4(1) / corner4(3));
    corner4_x = floor(corner4(2) / corner4(3));
    
    corner_y = [corner1_y, corner2_y, corner3_y, corner4_y]
    corner_x = [corner1_x, corner2_x, corner3_x, corner4_x]
    
    ymin = min(corner_y);
    ymax = max(corner_y);
    xmin = min(corner_y);
    xmax = max(corner_x);
    
    n_x_points = xmax - xmin;
    n_y_points = ymax - ymin;  
    
    % Generate a regular grid of (x,y) coordinates containing the output
    % image
    [x,y] = meshgrid(linspace(xmin,xmax,n_x_points), linspace(ymin,ymax,n_y_points));
   
    % Allocate an array of zeros with the same dimensions as your grid in
    % which you will store the output intensities
    [m, ~] = size(x(:));
    X = [x(:) y(:) ones(m,1)];

    % Apply inverse of H to the grid of coordinates to find the
    % corresponding locations in I
    Z = H \ X';
    Z = Z';
    
    Z1 = Z(:,1) ./ Z(:,3);
    Z2 = Z(:,2) ./ Z(:,3);
    Xh = [Z1 Z2];
    
    % Sample I at these points, using interp2 in the three channels
    if nChannels == 3
        Ih1 = interp2(double(I(:,:,1)), Xh(:,1), Xh(:,2), 'linear');
        Ih2 = interp2(double(I(:,:,2)), Xh(:,1), Xh(:,2), 'linear');
        Ih3 = interp2(double(I(:,:,3)), Xh(:,1), Xh(:,2), 'linear');
        Ih = [Ih1 Ih2 Ih3];
    end
    if nChannels == 1
        Ih = interp2(double(I(:,:,1)), Xh(:,1), Xh(:,2), 'linear');
    end

    % Reshape intensity vector with new height and width values
    Ih = reshape(Ih, [n_y_points, n_x_points, nChannels]);

    % Set intensities from outside boundaries to zero
    Ih(isnan(Ih)) = 0;
    
end

