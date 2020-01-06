function H = homography2d(x1,x2)
% Calculate homography matrix considering 4 points with correspondences
[Ncoords, Npoints] = size(x1);

if Npoints < 4
    error('Invalid no. of correspondences for homography estimation');
end

% Find the centroid of the correspondences
xcentroid1 = 0;
ycentroid1 = 0;
xcentroid2 = 0;
ycentroid2 = 0;

for i = 1:Npoints
    
    xcentroid1 = xcentroid1 + x1(1,i);
    ycentroid1 = ycentroid1 + x1(2,i);
    xcentroid2 = xcentroid2 + x2(1,i);
    ycentroid2 = ycentroid2 + x2(2,i);
    
end
xcentroid1 = xcentroid1/Npoints;
ycentroid1 = ycentroid1/Npoints;
xcentroid2 = xcentroid2/Npoints;
ycentroid2 = ycentroid2/Npoints;

% Generating the Translation Matrices
tr1 = [1, 0, -xcentroid1; 0, 1, -ycentroid1; 0, 0, 1];
tr2 = [1, 0, -xcentroid2; 0, 1, -ycentroid2; 0, 0, 1];

% Generating Translated Correspondences
x1_trans = zeros(Ncoords,Npoints);
x2_trans = zeros(Ncoords,Npoints);

for i = 1:Npoints
    
    x1_trans(:,i) = tr1*[x1(1,i); x1(2,i); 1];
    x2_trans(:,i) = tr2*[x2(1,i); x2(2,i); 1];
    
end
    
% Computing the average RMS distance of translated correspondences
rms1 = 0;
rms2 = 0;

for i = 1:Npoints 
    rms1 = rms1 + sqrt(x1_trans(1,i)^2 + x1_trans(2,i)^2);
    rms2 = rms2 + sqrt(x2_trans(1,i)^2 + x2_trans(2,i)^2);  
end

rms1 = rms1/Npoints;
rms2 = rms2/Npoints;

% Generating the scaling matrices
scale_factor1 = sqrt(2)/rms1;
scale_factor2 = sqrt(2)/rms2;
sc1 = [scale_factor1, 0, 0; 0, scale_factor1, 0; 0, 0, 1];
sc2 = [scale_factor2, 0, 0; 0, scale_factor2, 0; 0, 0, 1];

% Normalizing correspondences to average RMS distance sqrt(2)
x1_norm = zeros(Ncoords,Npoints);
x2_norm = zeros(Ncoords,Npoints);

for i = 1:Npoints
    
    x1_norm(:,i) = sc1*[x1_trans(1,i); x1_trans(2,i); 1];
    x2_norm(:,i) = sc2*[x2_trans(1,i); x2_trans(2,i); 1];
    
end
% Generating Transformation matrices and required inverse for
% denormalization
T1 = sc1*tr1;
T2 = sc2*tr2;

% Generation of the A matrix for computation H
A = zeros(2*Npoints, 9);

for i = 1:Npoints
    A(i,:) = [-x1_norm(1,i) -x1_norm(2,i) -x1_norm(3,i) 0 0 0 x2_norm(1,i)*x1_norm(1,i) x2_norm(1,i)*x1_norm(2,i) x2_norm(1,i)];
    A(2*i,:) = [0 0 0 -x1_norm(1,i) -x1_norm(2,i) -x1_norm(3,i) x2_norm(2,i)*x1_norm(1,i) x2_norm(2,i)*x1_norm(2,i) x2_norm(2,i)];
end

% Solve linear equation Ah = 0 using the SVD
if ((sum(sum(isnan(A))) + (sum(sum(isinf(A))) ~= 0)))
    S = zeros(1,9);
else
    [~, ~, V] = svd(A, 0);
    S = V(:,9);
end

H = [S(1) S(2) S(3); S(4) S(5) S(6); S(7) S(8) S(9)];

% Denormalizing H matrix
H = (T1\H)*T2;

end

