%% Author: M.K,  Winter 2017
% SeamCarving Algorithm
%

clc
close all
% in_img = imread ('ocean.bmp');
% in_img = imread('edo.png');
 in_img = imread('van.bmp');
% img = imread('tree.jpg');
% in_img = imread('bridge.jpg');
% img = imread('waterfall.png');
% in_img = imread('tower.jpg');
figure; imshow(in_img), title(' Original Image ');
seamsNum = input('\n Enter the number of seamsNum to be carved/expanded: \n');
carv_exp = input('\n Enter "0" to CARVE or "1" to EXPAND the seamsNum: \n');
hv = input('\n Enter "0" for HORIZONTAL seamsNum: (Anything else for Vertical) \n');
img = in_img;
if hv == 0
    img = imrotate(img,90);
end
% carve seamsNum
if carv_exp == 0
    im2 = double( img );
    for i = 1 : seamsNum
        energy = energyFunc( img );
        lowEnSeams = findseam( energy );
        [ fin_im , im2 ] = carveSeams( img , lowEnSeams , im2 );
        img = uint8( fin_im );
    end
    
 sn = [];   
% expand image
elseif carv_exp == 1
    if seamsNum > 50
        fin_im = img;
        temp = seamsNum ;
        while temp > 50
            sn1 = 50 ;
            temp = temp - sn1 ;
            energy = energyFunc( img );
            lowEnSeams  = findseam ( energy );
            [ fin_im , im2 ] = expandSeams( fin_im , lowEnSeams , sn1);
        end
    else
        energy = energyFunc( img );
        [ lowEnSeams ] = findseam ( energy );
        [ fin_im , im2 ] = expandSeams ( img, lowEnSeams , seamsNum );
    end
end
if hv == 0
    fin_im = imrotate( fin_im , 270 );
    im2 = imrotate( im2 , 270 );
end

if seamsNum > 50 && carv_exp == 1
    figure; imshow(uint8(fin_im)),title('Modified Image'); 
    imwrite(uint8(fin_im), 'Modified_image.bmp');
else
    figure; imshow(uint8(fin_im)),title('Modified Image'); 
    imwrite(uint8(fin_im), 'Modified_image.bmp');
    figure; imshow(uint8(im2)), title('Seams');
    imwrite(uint8(im2), 'seamsNum.bmp');
end
size(fin_im)

%% find the seamsNum with lowest energy
function  lowEnSeams  = findseam(energy) 
matrx = double(energy);
% add inf value to the left and right columns
%matrx = concat(matrx, [1 inf]);
matrx = padarray( matrx , [0,1], inf ); 
[m , n , p] = size(matrx);
for i = 2 : m 
    for j = 2 : ( n - 1)
        % find the smallest neighbor in the top 3 above each pixel
        neighbors = [matrx(i - 1, j - 1) matrx(i - 1, j) matrx(i - 1, j + 1)] ;
        % cumulative min energy for all possible connected seams
        matrx(i, j) = matrx(i, j) + min( neighbors );
    end
end
% find the min element in the last row
[val, inx] = min( matrx ( i , :) );     
seamEnergy = double( val );
lowEnSeams = ones( size( energy ) );
% find seam mask from the above matrix 'matrx' with lowest energy
    for i = m : -1 : 2
        % each pixel which in the found seams, its corresponding value in
        % the mask will be 0, otherwise 1
        lowEnSeams( i , inx - 1 ) = 0;
        neighbour = [ matrx(i-1,inx-1), matrx(i-1,inx), matrx(i-1,inx+1)];
        [~ , ind2 ] = min( neighbour );
        % seamEnergy = seamEnergy + double(minval);
        inx = inx + (ind2 - 2);
    end
end

%% function for carving the seamsNum with lowest energy
function [carvedIm,im2] = carveSeams(img, lowEnSeams,im2) % to remove the seam from the findseam
carvedIm = zeros(size(img, 1), size(img, 2) -1, size(img, 3)); % the actual new image size
drawSeams = zeros(size(im2, 1), size(im2, 2) , size(im2, 3)); % the image with seam shown as red
for i = 1 : size(lowEnSeams, 1)
    for j = find(lowEnSeams(i, :) == 0 )
        carvedIm(i, :, 1) = [img(i, 1:j-1, 1), img(i, j+1:end, 1) ];
        carvedIm(i, :, 2) = [img(i, 1:j-1, 2), img(i, j+1:end, 2) ]; 
        carvedIm(i, :, 3) = [img(i, 1:j-1, 3), img(i, j+1:end, 3) ];
        drawSeams(i,:,1) = [ im2(i,1:j-1,1) , 0 , im2(i,j+1:end,1) ];
        drawSeams(i,:,2) = [ im2(i,1:j-1,2) , 0 , im2(i,j+1:end,2) ];
        drawSeams(i,:,3) = [ im2(i,1:j-1,3) , 255 , im2(i,j+1:end,3) ];
    end
end
im2 = drawSeams;
end

%% the Energy Function
function res = energyFunc( image )
%     image = image(:,:,1);
%     filt = [ -1, 0, 1 ] ;
%     res = abs( imfilter( image, filt, 'replicate')) + abs( imfilter( image,filt' , 'replicate'));
In = im2double(image);
    grad1 = abs(gradient(In(:,:,1)));
    GImg2 = abs(gradient(In(:,:,2)));
    GImg3 = abs(gradient(In(:,:,3)));
    res = grad1 + GImg2 + GImg3;
end

%% function to expand the seamsNum
function [expandedIm, im2] = expandSeams(img, lowEnSeams , seamsNum)
img = double(img);
im2 = double(img);
% interpol = @(img, i, j, k) (img(i, j-1, k) + img(i, j+1, k))/2;
while seamsNum > 0
    expandedIm = zeros(size(img, 1), size(img, 2) +1 , size(img, 3));
    drawSeams = zeros(size(img, 1), size(img, 2) + 1, size(img, 3));
    for i = 1 : size(lowEnSeams, 1)
  % value of 0 means it is a seam
        for j = find(lowEnSeams(i, :) == 0)
            if j == size(lowEnSeams, 2)
                expandedIm(i, :, 1) = [img(i, 1:j, 1), img(i, j, 1), img(i, j+1:end, 1)];
                expandedIm(i, :, 2) = [img(i, 1:j, 2), img(i, j, 2), img(i, j+1:end, 2)];
                expandedIm(i, :, 3) = [img(i, 1:j, 3), img(i, j, 3), img(i, j+1:end, 3)];
                drawSeams(i,:,1) = [im2(i,1:j-1,1) , [0 0] , im2(i,j+1:end,1)];
                drawSeams(i,:,2) = [im2(i,1:j-1,2) , [0 0] , im2(i,j+1:end,2)];
                drawSeams(i,:,3) = [im2(i,1:j-1,3) , [255 255] , im2(i,j+1:end,3)];
            else
                expandedIm(i, :, 1) = [img(i, 1:j, 1), round(img(i, j, 1)./2 + img(i, j+1, 1)./2), img(i, j+1:end, 1)];
                expandedIm(i, :, 2) = [img(i, 1:j, 2), round(img(i, j, 2)./2 + img(i, j+1, 2)./2), img(i, j+1:end, 2)];
                expandedIm(i, :, 3) = [img(i, 1:j, 3), round(img(i, j, 3)./2 + img(i, j+1, 3)./2), img(i, j+1:end, 3)];       
                drawSeams(i,:,1) = [ im2( i , 1:j-1 , 1 ) , [0,0] , im2( i , j+1:end , 1 ) ];
                drawSeams(i,:,2) = [ im2( i , 1:j-1 , 2 ) , [0,0], im2( i , j+1:end , 2 ) ];
                drawSeams(i,:,3) = [ im2( i , 1:j-1 , 3 ) , [255,255] , im2( i , j+1:end , 3 ) ];
            end
        end
    end
    img = uint8( expandedIm );
    im2 = drawSeams ;
    energy = energyFunc( im2 );
    lowEnSeams = findseam( energy );
    seamsNum = seamsNum - 1 ;
end
end
