function [H, idx_inliers] = ransac_homography_adaptive_loop(x1, x2, th, max_it)

[Ncoords, Npoints] = size(x1);

% ransac
it = 0;
best_inliers = [];

while it < max_it
    
    points = randomsample(Npoints, 4);
    H = homography2d(x1(:,points), x2(:,points)); % ToDo: you have to create this function
    inliers = compute_inliers(H, x1, x2, th);
    
    % test if it is the best model so far
    if length(inliers) > length(best_inliers)
        best_inliers = inliers;
        [~,Ninliers] = size(best_inliers);
        Ninliers
    end    
    
    % update estimate of max_it (the number of trials) to ensure we pick, 
    % with probability p, an initial data set with no outliers
    fracinliers =  length(inliers)/Npoints;
    pNoOutliers = 1 -  fracinliers^4;
    pNoOutliers = max(eps, pNoOutliers);  % avoid division by -Inf
    pNoOutliers = min(1-eps, pNoOutliers);% avoid division by 0
    p=0.99;
    max_it = log(1-p)/log(pNoOutliers);
    
    it = it + 1;
end

% compute H from all the inliers
H = homography2d(x1(:,best_inliers), x2(:,best_inliers));
%H = goldStd2D(x1(:,best_inliers),x2(:,best_inliers));

idx_inliers = best_inliers;


function idx_inliers = compute_inliers(H, x1, x2, th)
    % Check that H is invertible
    if abs(log(cond(H))) > 15
        idx_inliers = [];
        return
    end
    
    % compute the symmetric geometric error % ToDo
    [Ncoords, Npoints] = size(x1);
    d2 = zeros(1,Npoints);
    
    for i = 1:Npoints
        
        x1_correspondence = H*[x1(1,i); x1(2,i); 1];
        x1_correspondence = x1_correspondence/x1_correspondence(3);
        
        x2_correspondence = H\[x2(1,i); x2(2,i); 1];
        x2_correspondence = x2_correspondence/x2_correspondence(3);
        
        error1 = sum(abs([x1(1,i); x1(2,i); 1] - x1_correspondence));
        error2 = sum(abs([x2(1,i); x2(2,i); 1] - x2_correspondence));
        
        d2(1,i) = error1 + error2;
    end
    
       idx_inliers = find(d2 < th.^2);

function xn = normalise(x)    
    xn = x ./ repmat(x(end,:), size(x,1), 1);

    
function item = randomsample(npts, n)
	a = [1:npts]; 
    item = zeros(1,n);    
    for i = 1:n
	  % Generate random value in the appropriate range 
	  r = ceil((npts-i+1).*rand);
	  item(i) = a(r);       % Select the rth element from the list
	  a(r)    = a(end-i+1); % Overwrite selected element
    end                       % ... and repeat