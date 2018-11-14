%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A(1) Fill in the missing values, given the partial values of the parameters 
%of the left and right cameras
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Image L
 % You are given
 % The projection matrix
    M_L= [1100.504780,0,331.023000,0; 0, 1097.763735,259.386377, 0; 0, 0,1,0]
    
 %Rotation: 
     R_L=[ 1  0 0 
         0 1 0
         0 0 1]
%Focal length: 
    f_L=1.0
    
% Compute COP1:
  T_L=

%Compute the image center: 
     ox_L = 
     oy_L =  
     
% compute the scale factor:
     Sx_L =  
     Sy_L = 
     
 % Compute intrinsic projection matrices: 

    Mint_L=    
     
% Image R: 
%You are given the intrinsic parameters: 
    %Image center: 
    ox_R = 320.798101, 
    oy_R = 236.431326,

    % Scale factor: 
    Sx_R = 1095.671499, 
    Sy_R = 1094.559584 
   
    % Focal length: 
    f_R=1
     
    % Translation w.r.t. the world origin
    T_R=-[178.2218
          18.8171
          -13.7744]
    % Rotation
    R_R =[
         0.9891    0.0602   -0.1346
        -0.0590    0.9982    0.0134
         0.1351   -0.0053    0.9908]

% Compute intrinsic projection matrices: Mint1 and Mint2

   Mint_R=
      
% Compute the projection matrix
    M_R=

% Compute 
    COP_L = 
    COP_R = 
    
% What is the distance between COP_L and COP_R?
D=
    
    
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
    f1=gcf;
    
    figure;
    imshow(im_R,[])
    hold on;
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
    plot(q_L(1),q_R(2),'*r');
    
    figure(f4);
    plot(q_R(1),q_R(2),'*r');
    
  
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A(5-7)compute the fundamnral matrix F, the epipoles e_L   % 
% and e_R,                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   
 
     e_L=
     e_R=
 
     F=
     
 % Please normalize F by F(3,3).
 
   
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % A(8) Draw epipolar lines                           
 % You have to write the function
 % [f1,f2]=draw_epipolar_lines(im_L,im_R,F,p,f1,f2) (f1 and f2 are the figures on
 % which im1 and im2 are already dispalyed, the value of f1 and f2 are 0 if the images
 % are note displayed yet. In this case, the f1 and f2 are the output of
 % the function for later use of the figures
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 p=[283,348]';% the selected point
 
 draw_epipolar_lines(im_L,im_R,F,p,f1,f2)

 figure(f1);
 
 % Choose points from image 1 (look at help getpts)
 figure(f1);
 [Px,Py]=getpts
 
% Display the  set of pipolar lines which corresponds to the chosen points

for i=1:length(Px)
    draw_epipolar_lines(im_L,im_R,F,[Px(i),Py(i)]',f1,f2)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % A(9) Featute detection, matching, and remove outliers
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
 Points_L=  detectSURFFeatures(im_L);
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
 title('100 Strongest Feature Points from  left image');
 hold on;
 plot(selectStrongest(Points_R, num_points));
   
 
 % Consider feature with index ind_L in the left image.
 % Write a function 'match_best()' that finds the most similar feature to Features_L(ind_L)
 % in the right image. The list of fetaures in the right image is given by
 % Features_R.
 % Test it for ind_L=10.
 index_R=match_best(Features_L, Features_R, indel_L)
 
 % Use matlab function for matching all the features (instead of applying
 % the oned you wrote)
 
 Pairs = matchFeatures(Features_L, Features_R);
 
 showMatchedFeatures(im_L, im_R,Points_L(Pairs(:,1),:),Points_R(Pairs(:,2),:),'montage');
 
 % Write a function that remove incorrect matches using epipolar geometry. 
 %  That is, any pairs that is not consistent with the epipolar geometry should be removed.
 % Note that we cannot test whether  p^T Fq=0 since it will never be exactly 0 due to noise.
% Hence, you should compute how far the point q is from the epipolar line given by p^T F.
% Use the function  'sampsonDistance()' for computing the distance between a point and a
 % line
 
 Pairs_clean=remove_incorrect_matches(Pairs,Points_L,Points_R,F,th)
 
 % display the matching features as before with only the correct mathced features .
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calibrated Stereo - Compute Matching
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ADD YOUR PART HERE



 
  