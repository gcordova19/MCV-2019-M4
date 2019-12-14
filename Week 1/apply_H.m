function Ih = apply_H(I, H) 
    %%%%%%%%%%%%%%%%%%%%
    % help.pdf to understand the process
    % STILL NOT WORKING
    %%%%%%%%%%%%%%%%%%%%%
    
    % Apply H to corners of I to determine the horizontal and vertical 
    % limits of the output image

    [width, height] = size(I);
    
    min = H \ [1; 1; 1];
    xmin = min(1) / min(3);
    ymin = min(2) / min(3);
    
    max = H \ [width; height; 1];
    xmax = max(1) / max(3);
    ymax = max(2) / max(3);
    
    % Generate a regular grid of (x,y) coordinates containing the output
    % image
    [x,y] = meshgrid(linspace(xmin,xmax,width), linspace(ymin,ymax,height));

    % Allocate an array of zeros with the same dimensions as your grid in
    % which you will store the output intensities
    [m, n] = size(x(:));
    X = [x(:) y(:) ones(m,1)];

    % Apply inverse of H to the grid of coordinates to find the
    % corresponding locations in I
    Z = H \ X';
    Z = Z';
    
    Z1 = Z(:,1) ./ Z(:,3);
    Z2 = Z(:,2) ./ Z(:,3);
    Zh = [Z1 Z2];
    
    % ToDo compute non-homogeneous coordinates
    % ToDo save result in Xh Nx2 matrix
   
    % Sample I at these points, using interp2
    Ih = interp2(I, Xh(:,1), Xh(:,2), 'linear');
end

