function Ih = apply_H(I, H) 
    % Apply H to corners of I to determine the horizontal and vertical 
    % limits of the output image
    [height, width, nChannels] = size(I);
    corners = [0, 0, size(I,2), size(I,2);
               0, size(I,1), 0, size(I,1)];

    dime = size(corners,2);
    c3d = [corners;  ones(1,dime)];
    h2d = H * c3d;
    c2d = h2d(1:2,:)./ [h2d(3,:)' h2d(3,:)']';
    
    corners_x = c2d(1:2,:);
    
    minx = floor(min(corners_x(1,:)));
    maxx = ceil(max(corners_x(1,:)));
    miny = floor(min(corners_x(2,:)));
    maxy = ceil(max(corners_x(2,:)));
    
    n_x_points = maxx - minx;
    n_y_points = maxy - miny;  
 
    % Generate a regular grid of (x,y) coordinates containing the output
    % image
    [x,y] = meshgrid(minx:maxx-1,miny:maxy-1);
   
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
    %Ih(isnan(Ih)) = 0;
    
end
