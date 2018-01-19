function practical2b

%The goal of this part of the practical is to take a real image containing
%a planar black square and figure out the transformation between the square
%and the camera.  We will then draw a wire-frame cube with it's base
%corners at the corner of the square.  You should use this
%template for your code and fill in the missing sections marked "TO DO"

%load in image 
im = imread('test104.jpg');

%define points on image
xImCart = [  140.3464  212.1129  346.3065  298.1344   247.9962;...
             308.9825  236.7646  255.4416  340.7335   281.5895];
         
%define 3D points of plane
XCart = [-50 -50  50  50 0 ;...
          50 -50 -50  50 0;...
           0   0   0   0 0];

%We assume that the intrinsic camera matrix K is known and has values
K = [640  0    320;...
     0    640  240;
     0    0    1];

%draw image and 2d points
figure; set(gcf,'Color',[1 1 1]);

       
%TO DO Use your routine to calculate TEst, the extrinsic matrix relating the
%plane position to the camera position.
TEst=estimatePlanePose(xImCart,XCart,K)



%define 3D points of plane
XWireFrameCart = [-50 -50  50  50 -50 -50  50  50;...
                   50 -50 -50  50  50 -50 -50  50;...
                    0   0   0   0 -100 -100 -100 -100];

%TO DO Draw a wire frame cube, by projecting the vertices of a 3D cube
%through the projective camera and drawing lines betweeen the resulting 2d image
%points
XImWireFrameCart=projectiveCamera(K,TEst,XWireFrameCart)
imshow(im); axis off; axis image; hold on;
%plot(cube_point_xImCart(1,:),cube_point_xImCart(2,:),'r.','MarkerSize',10);
num = size(XImWireFrameCart,2);
for cPoint = 1 : num
    if mod(cPoint, num / 2) ~= 0
        plot(XImWireFrameCart(1,cPoint:cPoint+1), XImWireFrameCart(2,cPoint:cPoint+1), 'r-', 'LineWidth',2);
    elseif cPoint == num / 2
        plot([XImWireFrameCart(1,cPoint), XImWireFrameCart(1,1)], [XImWireFrameCart(2,cPoint), XImWireFrameCart(2,1)], 'c-o', 'LineWidth',2);
    elseif cPoint == num
        plot([XImWireFrameCart(1,cPoint), XImWireFrameCart(1,1 + num/2)], [XImWireFrameCart(2,cPoint), XImWireFrameCart(2,1 + num/2)], 'g-', 'LineWidth',2);
    end
end
for cPoint = 1 : num / 2
    plot([XImWireFrameCart(1,cPoint), XImWireFrameCart(1,cPoint+4)],[XImWireFrameCart(2,cPoint), XImWireFrameCart(2,cPoint+4)],'b-', 'LineWidth',2);
end
%QUESTIONS TO THINK ABOUT...

%Do the results look realistic?
%If not, then what factors do you think might be causing this?


function xImCart = projectiveCamera(K,T,XCart);

%replace this
xImCart = [];

%TO DO convert Cartesian 3d points XCart to homogeneous coordinates XHom
XHom=[XCart; ones(1, size(XCart,2))]
%TO DO apply extrinsic matrix to XHom to move to frame of reference of
%camera
xCamHom=T*XHom;
%TO DO project points into normalized camera coordinates xCamHom by (achieved by
%removing fourth row)
xCamHom(4,:)=[];
%TO DO move points to image coordinates xImHom by applying intrinsic matrix
xImHom=K*xCamHom;
%TO DO convert points back to Cartesian coordinates xImCart
xImCart=xImHom(1:2,:)./repmat(xImHom(3,:),2,1);


function T = estimatePlanePose(xImCart,XCart,K)

%replace this
T = [];

%TO DO Convert Cartesian image points xImCart to homogeneous representation
%xImHom
xImHom=[xImCart; ones(1,size(xImCart,2))];
%TO DO Convert image co-ordinates xImHom to normalized camera coordinates
%xCamHom
xCamHom=inv(K)*xImHom;

%TO DO Estimate homography H mapping homogeneous (x,y)
%coordinates of positions in real world to xCamHom.  Use the routine you wrote for
%Practical 1B.
H=calcBestHomography(XCart, xCamHom);
%%TO DO Estimate first two columns of rotation matrix R from the first two
%columns of H using the SVD
[U, S, V] = svd(H(:,1:2));
R_12 = U * [1 0; 0 1; 0 0]*V';

%TO DO Estimate the third column of the rotation matrix by taking the cross
%product of the first two columns
R_3=cross(R_12(:,1) , R_12(:,2));

%TO DO Check that the determinant of the rotation matrix is positive - if
%not then multiply last column by -1.
R=[R_12 R_3];

if det(R)<=0
    %R(3,:)=(-1)*R(:,3);
    R_3=-1*R_3;
end
%TO DO Estimate the translation t by finding the appropriate scaling factor k
%and applying it to the third colulmn of H
k=sum(sum(R(:,1:2)./H(:,1:2)))/6;
Translation=k.*H(:,end);
%TO DO Check whether t_z is negative - if it is then multiply t by -1 and
%the first two columns of R by -1.
if Translation(end) < 0
    Translation = Translation * -1;
    R(:,1:2) = R(:,1:2)*-1;
end 
%assemble transformation into matrix form
T  = [R Translation;0 0 0 1];

function H = calcBestHomography(pts1Cart, pts2Cart)

    %should apply direct linear transform (DLT) algorithm to calculate best
    %homography that maps the points in pts1Cart to their corresonding matchin in 
    %pts2Cart

    %****TO DO ****: replace this
    %H = eye(3);

    %**** TO DO ****;
    %first turn points to homogeneous
    pts1_H = [pts1Cart; ones(1,size(pts1Cart,2))]
    pts2_H = [pts2Cart; ones(1,size(pts2Cart,2))]
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
        submatx_2(2*i-1,:)=pts1_H(i,:);
        submatx_3(2*i-1,:)=pts1_H(i,:)*(-pts2_H(i,2));
        submatx_3(2*i,:)=pts1_H(i,:)*(-pts2_H(i,1));
    end

    A=[submatx_1 submatx_2 submatx_3];
    %solve Ah = 0 by calling
    h = solveAXEqualsZero(A); %(you have to write this routine too - see below)

    %reshape h into the matrix H
    H_transpose = reshape(h,4,3);
    H = H_transpose';

function x = solveAXEqualsZero(A);
    [U,S,V] = svd(A);
    x = V(:,size(V,2));
