function practical1B

%the aim of the second part of practical 1 is to use the homography routine
%that you established in the first part of the practical.  We are going to
%make a panorama of several images that are related by a homography.  I
%provide 3 images (one of which is has a large surrounding region) and a
%matching set of points between these images.

%close all open figures
close all;

%load in the required data
load('PracticalData','im1','im2','im3','pts1','pts2','pts3','pts1b');
%im1 is center image with grey background
%im2 is left image 
%pts1 and pts2 are matching points between image1 and image2
%im3 is right image
%pts1b and pts3 are matching points between image 1 and image 3

%show images and points
figure; set(gcf,'Color',[1 1 1]);image(uint8(im1));axis off;hold on;axis image;
plot(pts1(1,:),pts1(2,:),'r.'); 

plot(pts1b(1,:),pts1b(2,:),'m.');

figure; set(gcf,'Color',[1 1 1]);image(uint8(im2));axis off;hold on;axis image;
plot(pts2(1,:),pts2(2,:),'r.'); 
figure; set(gcf,'Color',[1 1 1]);image(uint8(im3));axis off;hold on;axis image;
plot(pts3(1,:),pts3(2,:),'m.'); 

%****TO DO**** 
%calculate homography from pts1 to pts2
Homo_1to2=calcBestHomography(pts1, pts2);
Homo_1to3=calcBestHomography(pts1b, pts3);
%****TO DO**** 
%for every pixel in image 1
    %transform this pixel position with your homography to find where it 
    %is in the coordinates of image 2
    %if it the transformed position is within the boundary of image 2 then 
        %copy pixel colour from image 2 pixel to current position in image 1 
        %draw new image1 (use drawnow to force it to draw)
    %end
%end;

[x1, y1, ~]=size(im1);
[x2, y2, ~]=size(im2);
[x3, y3, ~]=size(im3);

for row=1:x1
    for col=1:y1
        pos_pts2=Homo_1to2*[col row 1]';
        lambda=pos_pts2(3);
        im2_r=int64(round(pos_pts2(2)/lambda));
        im2_c=int64(round(pos_pts2(1)/lambda));
        if im2_r<=x2 && im2_r>0 && im2_c<=y2 && im2_c>0
            im1(row,col, :)= im2(im2_r, im2_c, :);
        end
    end
end


figure; imshow(uint8(im1));
        
%****TO DO****
%repeat the above process mapping image 3 to image 1.
for row=1:x1
    for col=1:y1
        pos_pts3=Homo_1to3*[col row 1]';
        %lambda=pos_pts2(3);
        im3_r=int64(round(pos_pts3(2)/pos_pts3(3)));
        im3_c=int64(round(pos_pts3(1)/pos_pts3(3)));
        if im3_r<=x3 && im3_r>0 && im3_c<=y3 && im3_c>0
            im1(row,col, :)= im3(im3_r, im3_c, :);
        end
    end
end

function H = calcBestHomography(pts1Cart, pts2Cart)

    %should apply direct linear transform (DLT) algorithm to calculate best
    %homography that maps the points in pts1Cart to their corresonding matchin in 
    %pts2Cart

    %****TO DO ****: replace this
    %H = eye(3);

    %**** TO DO ****;
    %first turn points to homogeneous
    pts1_H = [pts1Cart; ones(1,size(pts1Cart,2))];
    pts2_H = [pts2Cart; ones(1,size(pts2Cart,2))];
    %then construct A matrix which should be (10 x 9) in size
    [row, col] = size(pts1_H);
    A = zeros(col*2,row*3);
    pts1_H = pts1_H';
    pts2_H = pts2_H';
    submatx_1 = zeros(col * 2,row);
    submatx_2 = zeros(col * 2,row);
    submatx_3 = zeros(col * 2,row);

    for i=1:col
        submatx_1(2*i,:)=pts1_H(i,:);
        submatx_2(2*i-1,:)=-pts1_H(i,:);
        submatx_3(2*i-1,:)=pts1_H(i,:)*(pts2_H(i,2));
        submatx_3(2*i,:)=pts1_H(i,:)*(-pts2_H(i,1));
    end

    A=[submatx_1 submatx_2 submatx_3];
    %solve Ah = 0 by calling
    h = solveAXEqualsZero(A); %(you have to write this routine too - see below)

    %reshape h into the matrix H
    H_transpose = reshape(h,3,3);
    H = H_transpose';





function x = solveAXEqualsZero(A);
    [U,S,V] = svd(A);
    x = V(:,size(V,2));


