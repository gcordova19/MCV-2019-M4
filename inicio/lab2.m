%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lab 2: Image mosaics

addpath('/MATLAB Drive/inicio/sift');

% 1. Compute image correspondences
if exist('/MATLAB Drive/inicio/tmp.jpg', 'file')==2
  delete('/MATLAB Drive/inicio/tmp.jpg');
end
% Open images

% imargb = imread('Data/llanes/llanes_a.jpg');
% imbrgb = imread('Data/llanes/llanes_b.jpg');
% imcrgb = imread('Data/llanes/llanes_c.jpg');

imargb = imread('Data/castle_int/0016_s.png');
imbrgb = imread('Data/castle_int/0015_s.png');
imcrgb = imread('Data/castle_int/0014_s.png');

% imargb = imread('Data/aerial/site13/frame00000.png');
% imbrgb = imread('Data/aerial/site13/frame00002.png');
% imcrgb = imread('Data/aerial/site13/frame00003.png');
% imargb = double(imread('Data/aerial/site22/frame_00001.tif'))/255;
% imbrgb = double(imread('Data/aerial/site22/frame_00018.tif'))/255;
% imcrgb = double(imread('Data/aerial/site22/frame_00030.tif'))/255; 
ima = sum(double(imargb), 3) / 3 / 255;
imb = sum(double(imbrgb), 3) / 3 / 255;
imc = sum(double(imcrgb), 3) / 3 / 255;

% imargb = double(imread('Data/aerial/site22/frame_00001.tif'))/255;
% imbrgb = double(imread('Data/aerial/site22/frame_00018.tif'))/255;
% imcrgb = double(imread('Data/aerial/site22/frame_00030.tif'))/255;
%  im1 = imargb;
%  im2 = imbrgb;
% imc = imcrgb;
% tmpim1 = ima;
% tmpim2 = imb;
tmpim1 = imargb;
tmpim2 = imbrgb;
tmpim3 = imcrgb;


%% Compute SIFT keypoints
[points_a, desc_a] = sift(ima, 'Threshold', 0.01);
[points_b, desc_b] = sift(imb, 'Threshold', 0.01);
[points_c, desc_c] = sift(imc, 'Threshold', 0.01);

% figure;
% imshow(imargb);%image(imargb)
% hold on;
% plot(points_a(1,:), points_a(2,:),'+y');
% figure;
% imshow(imbrgb);%image(imbrgb);
% hold on;
% plot(points_b(1,:), points_b(2,:),'+y');
% figure;
% imshow(imcrgb);%image(imcrgb);
% hold on;
% plot(points_c(1,:), points_c(2,:),'+y');

%% Match SIFT keypoints 

% between a and b
matches_ab = siftmatch(desc_a, desc_b,5.5);
% figure;
% plotmatches(ima, imb, points_a(1:2,:), points_b(1:2,:), matches_ab, 'Stacking', 'v');

% between b and c
matches_bc = siftmatch(desc_b, desc_c);
% figure;
% plotmatches(imb, imc, points_b(1:2,:), points_c(1:2,:), matches_bc, 'Stacking', 'v');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Compute the homography (DLT algorithm) between image pairs

%% Compute homography (normalized DLT) between a and b, play with the homography
% th = 12;
% xab_a = [points_a(1:2, matches_ab(1,:)); ones(1, length(matches_ab))];
% xab_b = [points_b(1:2, matches_ab(2,:)); ones(1, length(matches_ab))];
% [Hab, inliers_ab] = ransac_homography_adaptive_loop(xab_a, xab_b, th, 1000); % ToDo: complete this function
% 
% figure;
% plotmatches(ima, imb, points_a(1:2,:), points_b(1:2,:), ...
%     matches_ab(:,inliers_ab), 'Stacking', 'v');
% 
% vgg_gui_H(imargb, imbrgb, Hab);


%% Compute homography (normalized DLT) between b and c, play with the homography
% xbc_b = [points_b(1:2, matches_bc(1,:)); ones(1, length(matches_bc))];
% xbc_c = [points_c(1:2, matches_bc(2,:)); ones(1, length(matches_bc))];
% [Hbc, inliers_bc] = ransac_homography_adaptive_loop(xbc_b, xbc_c, th, 1000); 
% 
% figure;
% plotmatches(imb, imc, points_b(1:2,:), points_c(1:2,:), ...
%     matches_bc(:,inliers_bc), 'Stacking', 'v');
% 
% vgg_gui_H(imbrgb, imcrgb, Hbc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Build the mosaic

% corners = [-400 1600 -100 1200];
% iwb = apply_H_v2(imbrgb, [1 0 0; 0 1 0; 0 0 1] , corners);   % ToDo: complete the call to the function
% iwa = apply_H_v2(imargb, Hab, corners);    % ToDo: complete the call to the function
% iwc = apply_H_v2(imcrgb, inv(Hbc), corners);    % ToDo: complete the call to the function
% 
% figure;
% imshow(max(iwc, max(iwb, iwa)));%image(max(iwc, max(iwb, iwa)));axis off;
% title('Mosaic A-B-C');

% ToDo: compute the mosaic with castle_int images
% ToDo: compute the mosaic with aerial images set 13
% ToDo: compute the mosaic with aerial images set 22
% ToDo: comment the results in every of the four cases: hypothetise why it works or
%       does not work

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. Refine the homography with the Gold Standard algorithm
x = ones(3,size(matches_ab,2));
xp = ones(3,size(matches_ab,2));
x(1:2,:) = [points_a(1,matches_ab(1,:))-1;points_a(2,matches_ab(1,:))-1];
xp(1:2,:) = [points_b(1,matches_ab(2,:))-1;points_b(2,matches_ab(2,:))-1];
% % Homography ab
%Hini = dlt2D(x,xp);
Hini = homography2d(x,xp);
%Y_initial = gs_errfunction( P0, Xobs );
Y_initial =sum(gs_errfunction( x, xp, Hini));
err_initial = sum( sum( Y_initial.^2 ));
% Minimize the function using Levenberg-Marquardt
options = optimset('Algorithm', 'levenberg-marquardt');
H = lsqnonlin(@(t) gs_errfunction(x,xp,t),Hini, [], [], options);
% % we show the geometric error before and after the refinement
f=gs_errfunction(x,xp,H);
err_final = sum( sum( f.^2 ));
fprintf(1, 'Gold standard reproj error initial %f, final %f\n', err_initial, err_final);

finalImage = stitching2im(tmpim1,tmpim2,H);
figure;
imshow(finalImage)
imwrite(finalImage,'/MATLAB Drive/inicio/tmp.jpg');
% figure;
% imshow(finalImage);
% title('Mosaic A-B with key points. A with yellow and B with blue','FontSize', 14);
% hold on;
% plot(x(1,:), x(2,:),'+y');
% plot(xp(1,:), xp(2,:),'+c');
% imwrite(finalImage,'/MATLAB Drive/inicio/Mosaic_a_b_kp.jpg');
imdrgb = imread('/MATLAB Drive/inicio/tmp.jpg');
tmpim4 = imdrgb;
figure;
imshow(tmpim4);
imd = sum(double(imdrgb), 3) / 3 / 255;
[points_d, desc_d] = sift(imd, 'Threshold', 0.01);
matches_cd = siftmatch(desc_c, desc_d,5.5);
xd = ones(3,size(matches_cd,2));
xpd = ones(3,size(matches_cd,2));
xd(1:2,:) = [points_c(1,matches_cd(1,:))-1;points_c(2,matches_cd(1,:))-1];
xpd(1:2,:) = [points_d(1,matches_cd(2,:))-1;points_d(2,matches_cd(2,:))-1];
% % Homography ab
Hinid = homography2d(xd,xpd);
%Hinid = dlt2D(xd,xpd);
Y_initiald =sum(gs_errfunction( xd, xpd, Hinid));
err_initiald = sum( sum( Y_initiald.^2 ));
% Minimize the function using Levenberg-Marquardt
Hd = lsqnonlin(@(t) gs_errfunction(xd,xpd,t),Hinid, [], [], options);
% % we show the geometric error before and after the refinement
fd=gs_errfunction(xd,xpd,Hd);
err_finald = sum( sum( fd.^2 ));
fprintf(1, 'Gold standard reproj error initial %f, final %f\n', err_initiald, err_finald);

finalImaged = stitching2im(tmpim4,tmpim3,Hd);
figure;
imshow(finalImaged);
imwrite(finalImaged,'/MATLAB Drive/inicio/tmp2.jpg');
% figure;
% imshow(finalImaged);
% hold on;
% plot(xd(1,:), xd(2,:),'+y');
% plot(xpd(1,:), xpd(2,:),'+c');
% title('Mosaic A-B-C with key points. A and B Stitched with yellow and C with blue','FontSize', 14);
% imwrite(finalImage,'/MATLAB Drive/inicio/Mosaic_a_b_c_kp.jpg');