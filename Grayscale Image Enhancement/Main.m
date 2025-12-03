clc;
clear;
close all;

addpath Set12;
i = imread('sample.png');
i_gt = imread('sample_gt.png');
I = i(:,:,1);
[counts, ~] =  imhist(I);
counts = counts';
[M,N] = size(I);
number_of_pixels = M*N;

% the grayscale histogram distribution of the original image
h = counts;
% the ideal uniform distribution histogram
u = number_of_pixels/256*ones(1,256);
% difference matrix
D = zeros(255,256);
for p=1:255
   for q=1:256
       if p==q
        D(p,q) = -1;
       elseif q==p+1
           D(p,q) = 1;
       else
           D(p,q) = 0;
       end
   end
end
CostFunction=@(x) Costfunc(x,h,u,D);% cost function

pop_size = 30; % population size
dim = 256; % problem dimensions
% search boundary
x_min = 0;
x_max = 255;
% max number of iterations
max_it = 500;

% ExpNum = 20;
% for exp=1:ExpNum
[~,best_pos_EBFGO] = BFGO(CostFunction,dim,pop_size,max_it,x_min,x_max);
best_pos_EBFGO = best_pos_EBFGO';

[~,best_pos_BFGO] = BFGO(CostFunction,dim,pop_size,max_it,x_min,x_max);
best_pos_BFGO = best_pos_BFGO';

[~,best_pos_PSO] = PSO_func(CostFunction,dim,pop_size,max_it,x_min,x_max);
best_pos_PSO = best_pos_PSO';

[~,best_pos_VPPSO] = VPPSO(CostFunction,dim,pop_size,max_it,x_min,x_max);
best_pos_VPPSO = best_pos_VPPSO';

[~,best_pos_SCA] = SCA_func(CostFunction,dim,pop_size,max_it,x_min,x_max);
best_pos_SCA = best_pos_SCA';

[~,best_pos_EDOLSCA] = EDOLSCA(CostFunction,dim,pop_size,max_it,x_min,x_max);
best_pos_EDOLSCA = best_pos_EDOLSCA';

%% calculate enhance result

lambda = 3;% control coefficient

c_he = zeros(256,1); 
c_he(1,1) = counts(1,1);
for n = 1:255
    c_he(n+1,1) = c_he(n,1) + counts(1,n+1);
end
for n = 0:255
    c_he(n+1,1) = c_he(n+1,1)/number_of_pixels;
end
cdf_he = round(255.*c_he + 0.5);

c_EBFGO = zeros(256,1); 
c_EBFGO(1,1) = best_pos_EBFGO(1,1);
for n = 1:255
    c_EBFGO(n+1,1) = c_EBFGO(n,1) + best_pos_EBFGO(n+1,1);
end
for n = 0:255
    c_EBFGO(n+1,1) = c_EBFGO(n+1,1)/number_of_pixels;
end
cdf_EBFGO = round(lambda*(255.*c_EBFGO + 0.5));

c_BFGO = zeros(256,1); 
c_BFGO(1,1) = best_pos_BFGO(1,1);
for n = 1:255
    c_BFGO(n+1,1) = c_BFGO(n,1) + best_pos_BFGO(n+1,1);
end
for n = 0:255
    c_BFGO(n+1,1) = c_BFGO(n+1,1)/number_of_pixels;
end
cdf_BFGO = round(lambda*(255.*c_BFGO + 0.5));

c_PSO = zeros(256,1); 
c_PSO(1,1) = best_pos_PSO(1,1);
for n = 1:255
    c_PSO(n+1,1) = c_PSO(n,1) + best_pos_PSO(n+1,1);
end
for n = 0:255
    c_PSO(n+1,1) = c_PSO(n+1,1)/number_of_pixels;
end
cdf_PSO = round(lambda*(255.*c_PSO + 0.5));

c_VPPSO = zeros(256,1); 
c_VPPSO(1,1) = best_pos_VPPSO(1,1);
for n = 1:255
    c_VPPSO(n+1,1) = c_VPPSO(n,1) + best_pos_VPPSO(n+1,1);
end
for n = 0:255
    c_VPPSO(n+1,1) = c_VPPSO(n+1,1)/number_of_pixels;
end
cdf_VPPSO = round(lambda*(255.*c_VPPSO + 0.5));

c_SCA = zeros(256,1); 
c_SCA(1,1) = best_pos_SCA(1,1);
for n = 1:255
    c_SCA(n+1,1) = c_SCA(n,1) + best_pos_SCA(n+1,1);
end
for n = 0:255
    c_SCA(n+1,1) = c_SCA(n+1,1)/number_of_pixels;
end
cdf_SCA = round(lambda*(255.*c_SCA + 0.5));

c_EDOLSCA = zeros(256,1); 
c_EDOLSCA(1,1) = best_pos_EDOLSCA(1,1);
for n = 1:255
    c_EDOLSCA(n+1,1) = c_EDOLSCA(n,1) + best_pos_EDOLSCA(n+1,1);
end
for n = 0:255
    c_EDOLSCA(n+1,1) = c_EDOLSCA(n+1,1)/number_of_pixels;
end
cdf_EDOLSCA = round(lambda*(255.*c_EDOLSCA + 0.5));

i_he = uint8(zeros(M,N,3));
i_EBFGO = uint8(zeros(M,N,3));
i_BFGO = uint8(zeros(M,N,3));
i_PSO = uint8(zeros(M,N,3));
i_VPPSO = uint8(zeros(M,N,3));
i_SCA = uint8(zeros(M,N,3));
i_EDOLSCA = uint8(zeros(M,N,3));
for l = 1:M
    for w = 1:N
         i_he(l,w,1) = cdf_he(I(l,w)+1,1);
         i_he(l,w,2) = i_he(l,w,1);
         i_he(l,w,3) = i_he(l,w,1);

         i_EBFGO(l,w,1) = cdf_EBFGO(I(l,w)+1,1);
         i_EBFGO(l,w,2) = i_EBFGO(l,w,1);
         i_EBFGO(l,w,3) = i_EBFGO(l,w,1);
         
         i_BFGO(l,w,1) = cdf_BFGO(I(l,w)+1,1);
         i_BFGO(l,w,2) = i_BFGO(l,w,1);
         i_BFGO(l,w,3) = i_BFGO(l,w,1);

         i_PSO(l,w,1) = cdf_PSO(I(l,w)+1,1);
         i_PSO(l,w,2) = i_PSO(l,w,1);
         i_PSO(l,w,3) = i_PSO(l,w,1);

         i_VPPSO(l,w,1) = cdf_VPPSO(I(l,w)+1,1);
         i_VPPSO(l,w,2) = i_VPPSO(l,w,1);
         i_VPPSO(l,w,3) = i_VPPSO(l,w,1);

         i_SCA(l,w,1) = cdf_SCA(I(l,w)+1,1);
         i_SCA(l,w,2) = i_SCA(l,w,1);
         i_SCA(l,w,3) = i_SCA(l,w,1);

         i_EDOLSCA(l,w,1) = cdf_EDOLSCA(I(l,w)+1,1);
         i_EDOLSCA(l,w,2) = i_EDOLSCA(l,w,1);
         i_EDOLSCA(l,w,3) = i_EDOLSCA(l,w,1);
    end
end
% display the fusion results
% figure, subplot(2,2,1); imhist(i_he); title('Histogram Equalization');
figure, imshow(i_gt),title('GT Image')
figure, imshow(i_he),title('Enhance Result of HE')
figure, imshow(i_EBFGO),title('Enhance Result of EBFGO')
figure, imshow(i_BFGO),title('Enhance Result of BFGO')
figure, imshow(i_VPPSO),title('Enhance Result of VPPSO')
figure, imshow(i_PSO),title('Enhance Result of PSO')
figure, imshow(i_EDOLSCA),title('Enhance Result of EDOLSCA')
figure, imshow(i_SCA),title('Enhance Result of SCA')
figure, imshow(i),title('Original Image')
% figure, hist(i_EBFGO(:,:,1),counts);% historam

% evaluating results

entropy_of_original_image = entropy(i);
entropy_of_he = entropy(i_he);
entropy_of_gt = entropy(i_gt);
entropy_of_EBFGO = entropy(i_EBFGO);
entropy_of_BFGO = entropy(i_BFGO);
entropy_of_PSO = entropy(i_PSO);
entropy_of_VPPSO = entropy(i_VPPSO);
entropy_of_SCA = entropy(i_SCA);
entropy_of_EDOLSCA = entropy(i_EDOLSCA);

% std2_org = std2(i);
% std2_he = std2(i_he);
% std2_EBFGO = std2(i_EBFGO);
% std2_BFGO = std2(i_BFGO);
% std2_PSO = std2(i_PSO);
% std2_VPPSO = std2(i_VPPSO);
% std2_SCA = std2(i_SCA);
% std2_EDOLSCA = std2(i_EDOLSCA);

psnr_org = psnr(i,i_gt);
psnr_he = psnr(i_he,i_gt);
psnr_EBFGO = psnr(i_EBFGO,i_gt);
psnr_BFGO = psnr(i_BFGO,i_gt);
psnr_PSO = psnr(i_PSO,i_gt);
psnr_VPPSO = psnr(i_VPPSO,i_gt);
psnr_SCA = psnr(i_SCA,i_gt);
psnr_EDOLSCA = psnr(i_EDOLSCA,i_gt);

ssim_org = ssim(i,i_gt);
ssim_he = ssim(i_he,i_gt);
ssim_EBFGO = ssim(i_EBFGO,i_gt);
ssim_BFGO = ssim(i_BFGO,i_gt);
ssim_PSO = ssim(i_PSO,i_gt);
ssim_VPPSO = ssim(i_VPPSO,i_gt);
ssim_SCA = ssim(i_SCA,i_gt);
ssim_EDOLSCA = ssim(i_EDOLSCA,i_gt);


% end
