%% SENSE Parallel Imaging
% Load data
clc;
clear all;
%Shepp-Logan Phantom - Object
p = phantom('Modified Shepp-Logan', 96);
input=load('data\profile0.mat');
imshow(p,[]);
raw=input.ffts;
%Calculating size of image
%Nx and Ny - dimensions of image
%Nc - Number of coils
[Nx,Ny,Nz,Nc] = size(raw);
show_grid(log(abs(raw)), [], jet);
img = (ifft2(raw));    
show_grid(abs(img),[], gray);
%%
%Computing square root of sum of squares of images

img_combined  = sqrt(sum(abs(img).^2,4));
show_img(abs(img_combined),[]);
S_0 = img./img_combined;
show_grid(abs(S_0),[],jet);
%% 
%Smoothing
kernel = ones(3)/3^2;
S_1 = zeros(Nx,Ny,Nz,Nc);
for i = 1:Nc
    S_1(:,:,1,i) = conv2(S_0(:,:,1,i),kernel,'same');
end
show_grid(abs(S_1),[],jet);
%Masking sensitivities
thresh = 0.0001*max(abs(img_combined(:)));
mask = abs(img_combined) > thresh;
%%
S_2 = S_1.*mask;
show_grid(abs(S_2),[],jet);
%% Use sensivities 

img_combined_opt = sum(img.*conj(S_2),4)./sum(S_2.*conj(S_2),4);
show_img(abs(img_combined_opt),[],gray);
%% SENSE Parallel Imaging

raw_R2 = raw;
raw_R2(1:2:end,:,:,:) = 0;
img_R2 = (ifft2(raw_R2));
show_grid(abs(img_R2),[],gray)
%% 

% initialise output image
img_R2_SENSE = zeros(Nx,Ny);
% loop over the top-half of the image
for x = 1:Nx/2
    % loop over the entire left-right extent
    for y = 1:Ny
        % pick out the sub-problem sensitivities
        S_R2 = transpose(reshape(S_2([x x+Nx/2],y,1,:),2,[]));
        % solve the sub-problem in the least-squares sense
        img_R2_SENSE([x x+Nx/2],y) = pinv(S_R2)*reshape(img_R2(x,y,1,:),[],1);
    end
end
% plot the result
show_img(abs(img_R2_SENSE),[],gray);

% Try solving the |R=3| problem:
%% 

raw_R3 = raw;
raw_R3(1:3:end,:,:,:) = 0;
raw_R3(2:3:end,:,:,:) = 0;
img_R3 = (ifft2(raw_R3));
show_grid(abs(img_R3),[],gray);
img_R3_SENSE = SENSE(img_R3,S_2,3);
show_img(abs(img_R3_SENSE),[],gray);
%% 
% and the |R=4| subproblem:

raw_R4 = raw;
raw_R4(1:4:end,:,:,:) = 0;
raw_R4(2:4:end,:,:,:) = 0;
raw_R4(3:4:end,:,:,:) = 0;
img_R4 = (ifft2(raw_R4));
show_grid(abs(img_R4),[],gray)
% define your solution to the R=4 problem
img_R4_SENSE = SENSE(img_R4,S_2,4);
show_img(abs(img_R4_SENSE),[],gray);
%% 
%% The Limits of Parallel Imaging
% Now that we have our SENSE implementation working, lets consider the limits 
% of parallel imaging. 
% 
% Try pushing the under-sampling factor to |R=6:|

% R=6 SENSE reconstruction
raw_R6 = raw;
idx_R6 = setdiff(1:Nx,6:6:Nx);
raw_R6(idx_R6,:,:,:) = 0;
img_R6 = ifft2(raw_R6);
img_R6_SENSE = SENSE(img_R6,S_2,6);
show_img(abs(img_R6_SENSE),[],gray);
%% 
% Or |R=8:|

% R=8 SENSE reconstruction
raw_R8 = raw;
idx_R8 = setdiff(1:Nx,8:8:Nx);
raw_R8(idx_R8,:,:,:) = 0;
img_R8 = ifft2(raw_R8);
img_R8_SENSE = SENSE(img_R8,S_2,8);
show_img(abs(img_R8_SENSE),[],gray);
%% 
% Along the same lines, what happens if you select a subset of the coil 
% information? Try a reconstruction at |R=4| using only |Nc=8| of the coil sensitivities 
% and corresponding images:

% Nc=8 SENSE reconstruction @ R=4
img_R4_Nc8 = SENSE(img_R4(:,:,:,1:8),S_2(:,:,:,1:8),4);
show_img(abs(img_R4_Nc8),[],gray);
%% 
% Or |Nc=4:|

% Nc=4 SENSE reconstruction @ R=4
img_R4_Nc4 = SENSE(img_R4(:,:,:,1:4),S_2(:,:,:,1:4),4);
show_img(abs(img_R4_Nc4),[],gray);
%% 
% * What happens to the reconstructions as |Nc| is close to or equal to |R? 
% |
% * How might you select an _optimal_ subset of |Nc=4| coils? What makes it 
% _optimal_?
% * What happens if you try to perform an |R=4| reconstruction with |Nc<R|?

% Optimal Nc=4 SENSE reconstruction @ R=4 using SVD-based coil compression
[~,~,v]=svd(reshape(S_2,[],16),0);
S_Nc4 = reshape(reshape(S_2,[],16)*v(:,1:4),[Nx,Ny,Nz,4]);
img_R4_Nc4 = reshape(reshape(img_R4,[],16)*v(:,1:4),[Nx,Ny,Nz,4]);
img_R4_Nc4opt = SENSE(img_R4_Nc4,S_Nc4,4);
show_img(abs(img_R4_Nc4opt),[0 16],gray);
%% 
% Noise SENSE reconstructions
outN4 = zeros(Nx,Ny,100);
outN1 = zeros(Nx,Ny,100);
for n = 1:100
    noise = (randn(Nx,Ny,Nz,Nc) + 1j*randn(Nx,Ny,Nz,Nc));
    outN4(:,:,n) = SENSE(noise,S_2,4);
    outN1(:,:,n) = SENSE(noise,S_2,1);
end
outN1(outN1==0)=inf; % prevents any division-by-zero
show_img(std(outN4,[],3)./(std(outN1,[],3)),[0 3],jet);
colorbar();
%% 
% * What do you observe happening to the noise of the output images as the under-sampling 
% factor increases? 
% 
% The spatial non-uniformity of the output noise is the so-called "geometry-factor", 
% which characterises noise-amplification due to the reconstruction. 
% 
% * Where is the noise amplification worst?
% * Can you relate this to some characteristic of the linear least squares sub-problems? 
% * How might you theoretically predict the amount of noise amplification you 
% expect to see?
% * How do you think this influences the choice of under-sampling factor?
%% Helper functions
%%
function out = SENSE(input,sens,R)
    [Nx,Ny,Nz,Nc] = size(input);
    out = zeros(Nx,Ny);
    % loop over the top-1/R of the image
    for x = 1:Nx/R
        x_idx = x:Nx/R:Nx;
        % loop over the entire left-right extent
        for y = 1:Ny
            % pick out the sub-problem sensitivities
            S = transpose(reshape(sens(x_idx,y,1,:),R,[]));
            % solve the sub-problem in the least-squares sense
            out(x_idx,y) = pinv(S)*reshape(input(x,y,1,:),[],1);
        end
    end
end


function show_grid(data, cscale, cmap)
    if nargin < 2
        cscale = [];
    end
    if nargin < 3
        cmap = gray;
    end
    figure();
    N = ndims(data);
    sz = size(data,N);
    for i = 1:sz
        subplot(sz / 4, 4, i);
        imshow(data(:, :, :, i), []);
    end
end

function show_img(data, cscale, cmap)
   if nargin < 2 || isempty(cscale)
       cscale = [-inf inf];
   end
   if nargin < 3
       cmap = gray;
   end
   figure();
   imagesc(data);
   axis equal
   colormap(cmap);
   caxis(cscale);
   plotH = gca;
   plotH.XTick = [];plotH.YTick = [];plotH.YColor = 'w';plotH.XColor = 'w';
end
function S = adaptive_est_sens(data)
    [Nx,Ny,Nz,Nc] = size(data);
    S = zeros(Nx,Ny,Nz,Nc);
    M = zeros(Nx,Ny,Nz);
    w = 5;
    for i = 1:Nx
        ii = max(i-w,1):min(i+w,Nx);
        for j = 1:Ny
            jj = max(j-w,1):min(j+w,Ny);
            for k = 1:Nz
                kk = max(k-w,1):min(k+w,Nz);
                kernel = reshape(data(ii,jj,kk,:),[],Nc);
                [V,D] = eigs(conj(kernel'*kernel),1);
                S(i,j,k,:) = V*exp(-1j*angle(V(1)));
                M(i,j,k) = sqrt(D);
            end
        end
    end
    S = S.*(M>0.1*max(abs(M(:))));
end
function out = fft2c(input)
    out = fftshift(fft(ifftshift(input,1),[],1),1);
    out = fftshift(fft(ifftshift(out,2),[],2),2);
end
function out = ifft2c(input)
    out = fftshift(ifft(ifftshift(input,1),[],1),1);
    out = fftshift(ifft(ifftshift(out,2),[],2),2);
end