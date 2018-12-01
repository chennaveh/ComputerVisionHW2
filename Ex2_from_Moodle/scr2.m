%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A(1) Fill in the missing values, given the partial values of the parameters 
%of the left and right cameras
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all;
% Image L
    % You are given
    % The projection matrix
    M_L= [1100.504780,0,331.023000,0; 0, 1097.763735,259.386377, 0; 0, 0,1,0];
    % M_L=[3.53553e+2, 3.39645e+2, 2.77744e+2, -1.44946e+6; -1.03528e+2, 2.33212e+1, 4.59607e+2, -6.32525e+5; 7.07107e-1 , -3.53553e-1, 6.12372e-1, -9.18559e+2];

    %Rotation: 
     R_L=[1 0 0 
          0 1 0
          0 0 1];
    %Focal length: 
    f_L=1.0;

    Null_vector = null(M_L);
    COP_L = Null_vector(1:3, :) / Null_vector(4);

    %Compute the image center:
    ox_L = M_L(1, 3);
    oy_L = M_L(2, 3);

    % compute the scale factor:
    Sx_L = M_L(1, 1);
    Sy_L = M_L(2, 2);

    % Compute intrinsic projection matrices: 
    Mint_L = M_L(:,1:3);
     
% Image R: 
%You are given the intrinsic parameters: 
    %Image center: 
    ox_R = 320.798101;
    oy_R = 236.431326;

    % Scale factor: 
    Sx_R = 1095.671499;
    Sy_R = 1094.559584;
   
    % Focal length: 
    f_R=1;
     
    % Translation w.r.t. the world origin
    T_R=-[178.2218
          18.8171
          -13.7744];
    % Rotation
    R_R =[0.9891    0.0602   -0.1346
         -0.0590    0.9982    0.0134
          0.1351   -0.0053    0.9908];

    % Compute intrinsic projection matrices: Mint1 and Mint2

    Mint_R = [Sx_R    0      ox_R
              0       Sy_R   oy_R
              0       0      1   ];
      
    % Compute the projection matrix
    Mext_R = eye(4);
    Mext_R(1:3,1:3) = R_R;
    Mext_R(1:3, 4) = -R_R*T_R; 
    aux = [eye(3), zeros([3 1])];

    M_R = Mint_R * aux * Mext_R;
    M_L= [1100.504780,0,331.023000,0; 0, 1097.763735,259.386377, 0; 0, 0,1,0];
    % Compute 
    % COP_L already calculated
    Null_vector = null(M_R);
    COP_R = Null_vector(1:3, :) / Null_vector(4);
    
    % What is the distance between COP_L and COP_R?
    D = norm(COP_R - COP_L);
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A(2) Projection                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
   % Write a finction p=proj(M,P) that recieves as input 
   % the 3D point P,  and a projection  matrix M, and output 
   % the  2D coordinates of the projected point.
    
    P=[-140,50,1200]';
    p_L=proj(M_L,P);
    p_R=proj(M_R,P);
    
    Q=[30,100,2000]'; 
    q_L=proj(M_L,Q);
    q_R=proj(M_R,Q);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read the images into im1 and im2,  and dispaly them.      %
% Using the command figure(f1) will redisplay the figure    %
% containing im1.                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
    im_L=imread('Left.tif');
    im_R=imread('Right.tif');

    figure;
    imshow(im_L,[])
    hold on;
    title('Left image with projection P');
    f1=gcf;
    
    figure;
    imshow(im_R,[])
    hold on;
    title('Rigth image with projection P');
    f2=gcf;

    % display the projection on the two images
    figure(f1);
    plot(p_L(1),p_L(2),'*r');
    
    figure(f2);
    plot(p_R(1),p_R(2),'*r');
    
    
    figure;
    imshow(im_L,[])
    hold on;
    f3=gcf;
    
    figure;
    imshow(im_R,[])
    hold on;
    f4=gcf;
    
    figure(f3);
    title('Left image with projection Q');
    plot(q_L(1),q_L(2),'*r');
    
    figure(f4);
    title('Right image with projection Q');
    plot(q_R(1),q_R(2),'*r');
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A(5-7)compute the fundamnral matrix F, the epipoles e_L   % 
% and e_R,                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Epipoles
    e_L = proj(M_L, COP_R);
    e_R = proj(M_R, COP_L);

% Fundamental matrix
    e_R_h = [e_R ; 1]; % homogeneous coordinates
    M_e_R = [0          -e_R_h(3)  e_R_h(2) 
             e_R_h(3)   0          -e_R_h(1) 
            -e_R_h(2)   e_R_h(1)   0 ]; % skewed symmetric matrix
        
    F = M_e_R * M_R * pinv(M_L);

% Please normalize F by F(3,3).
    F = F / F(3, 3);
    
    e_L_h = [e_L ; 1]; % homogeneous coordinates

% F-Epipoles consistency check
    % e_R_h' * F
    % F * e_L_h
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % A(8) Draw epipolar lines                           
 % You have to write the function
 % [f1,f2]=draw_epipolar_lines(im_L,im_R,F,p,f1,f2) (f1 and f2 are the figures on
 % which im1 and im2 are already dispalyed, the value of f1 and f2 are 0 if the images
 % are note displayed yet. In this case, the f1 and f2 are the output of
 % the function for later use of the figures
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    figure(f1);
 
% Choose points from image 1 (look at help getpts)
    figure(f1);
    [Px,Py]=getpts; %TODO - need save the points for the grader (see note at the end of the PDF)
    %pointsForGrader = [286.0000 160.0000];
    
% Display the  set of pipolar lines which corresponds to the chosen points

    for i=1:length(Px)
        draw_epipolar_lines(im_L,im_R,F,[Px(i),Py(i)]',f1,f2)
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % A(9) Featute detection, matching, and remove outliers
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
    Points_L = detectSURFFeatures(im_L);
    Points_R = detectSURFFeatures(im_R);
 
% Points_L & Points_R is an array of objects.
  
    [Features_L, Points_L] = extractFeatures(im_L,Points_L);
    [Features_R, Points_R] = extractFeatures(im_R,Points_R);
    num_points=100;

    figure;
    imshow(im_L);
    title('100 Strongest Feature Points from  left image');
    hold on;
    plot(selectStrongest(Points_L, num_points));

    figure;
    imshow(im_R);
    title('100 Strongest Feature Points from  left image');%TODO - shoud it be right image?
    hold on;
    plot(selectStrongest(Points_R, num_points));
   
 
 % Consider feature with index ind_L in the left image.
 % Write a function 'match_best()' that finds the most similar feature to Features_L(ind_L)
 % in the right image. The list of fetaures in the right image is given by
 % Features_R.
 % Test it for ind_L=10.
 
%     ind_L = 10;
%     ind_R = match_best(Features_L, Features_R, ind_L);
%     figure; imshow(im_L); title(['best match index = ' ind_L]);
%     hold on; plot(Points_L(ind_L));
%     
%     figure; imshow(im_R); title(['best match index = ' ind_R]);
%     hold on; plot(Points_R(ind_R));
    %%

 % Use matlab function for matching all the features (instead of applying
 % the oned you wrote)
 
    Pairs = matchFeatures(Features_L, Features_R);

    figure();
    subplot(2,1,1);
    showMatchedFeatures(im_L, im_R,Points_L(Pairs(:,1),:),Points_R(Pairs(:,2),:),'montage');
    title('Matching points');
 
% Write a function that remove incorrect matches using epipolar geometry. 
% That is, any pairs that is not consistent with the epipolar geometry should be removed.
% Note that we cannot test whether  p^T Fq=0 since it will never be exactly 0 due to noise.
% Hence, you should compute how far the point q is from the epipolar line given by p^T F.
% Use the function  'sampsonDistance()' for computing the distance between a point and a
% line
 
    th = 10;
    Pairs_clean = remove_incorrect_matches(Pairs,Points_L,Points_R,F,th);
    subplot(2,1,2);
    showMatchedFeatures(im_L, im_R,Points_L(Pairs_clean(:,1),:),Points_R(Pairs_clean(:,2),:),'montage');
% display the matching features as before with only the correct mathced features .
    title(['Matching points after incorrect points removed. th = ', num2str(th)]);
     

%% Part B

%clc; clear all; close all;

% Choose points from image left (look at help getpts)
figure;
imshow(im_L,[])
hold on;
title('Choose points on left image:')
f8=gcf;
[Px,Py]=getpts;
ps1 = [Px,Py];
%ps1 = [286.0000, 163.0000];

% Choose points from image right (look at help getpts)
figure;
imshow(im_R,[])
hold on;
title('Choose corresponded points on right image:')
f9=gcf;
[Px,Py]=getpts;
%ps2 = [314.0000 ,163.0000];
ps2 = [Px,Py];

P = stereo_list(ps1,ps2, M_L,M_R);

p_L = zeros(size(ps1,1),2);
p_R = zeros(size(ps2,1),2);
for i=1:size(P,1)
    p_L(i,:)=proj(M_L,P(i,:)');
    p_R(i,:)=proj(M_R,P(i,:)');
    
    %does p_L equal to ps1?
    errorLeft = abs(p_L(i,:)-ps1(i,:));
    %does p_R equal to ps2?
    errorRight = abs(p_R(i,:)-ps2(i,:));
end

figure(f8);
plot([p_L(:,1)],[p_L(:,2)],'*r');

figure(f9);
plot([p_R(:,1)],[p_R(:,2)],'*r');

%% Part C %%
    clear all;clc;
    im_L=imread('view5.tif');
    im_R=imread('view1.tif');
    
    T = 0.16; %distance between 2 cameras in meters
    
    D_out = disparityCalc(im_L,im_R,3,3,40,120);
    
    figure;
    imshow(D_out,[]);
    title('Disparity map view1 and view5');
    
    %C.e
    D2d = zeros(size(im_L,1),size(im_L,2),2);
    D2d(:,:,1) = D_out;
    d = imwarp(im_L,D2d);
    
    figure;
    imshow(d,[]);
    title('Back project view1 to view5');
    
    % C.f
    %Z = (f*T)/d -> T(distnace between cams), d(disparity) f=1
    Z = T*(ones(size(D_out)))./(D_out+100);
    figure;
    imshow(Z,[]);
    title('Depth map using disparity');
    
    % Qc.g - Repeat c-f 
    [Fx,Fy] = gradient(double(im_L));
    im_L_grad = (Fx.^2 + Fy.^2).^0.5;
    
    [Fx,Fy] = gradient(double(im_R));
    im_R_grad = (Fx.^2 + Fy.^2).^0.5;
        
    D_out = disparityCalc(im_L_grad,im_R_grad,3,3,40,120);
    figure;
    imshow(D_out,[]);
    title('Disparity map view1 and view5 with gradients');
    
     %C.e
    D2d = zeros(size(im_L,1),size(im_L,2),2);
    D2d(:,:,1) = D_out;
    d = imwarp(im_L,D2d);
        
    figure;
    imshow(d,[]);
    title('Back project view1 to view5 with gradients');
    
    % C.f
    %Z = (f*T)/d -> T(distnace between cams), d(disparity) f=1
    Z = T*(ones(size(D_out)))./(D_out+100);
    figure;
    imshow(Z,[]);
    title('Depth map with using gradients ');
    
    