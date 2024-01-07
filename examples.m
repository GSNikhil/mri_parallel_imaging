%   examples.m
%   mchiew@fmrib.ox.ac.uk

%   Load example data
input = matfile('data/data.mat');
truth = input.truth;
calib = input.calib;

%%{
input   =   matfile('data/profile0.mat');
truth   =   squeeze(input.ffts);
calib   =   squeeze(input.ffts);

truth = permute(truth, [3 2 1]);
calib = permute(calib, [3 2 1]);
%}

% Disp
figure;
title("Input FTs");
for i = 1:32
    subplot(8, 4, i);
    imshow(squeeze(log10(abs(truth(i, :, :)))), []);
end

% Disp
figure;
title("Input IFFTs");
for i = 1:32
    subplot(8, 4, i);
    a = squeeze(truth(i, :, :));
    a = (ifft2(a));
    imshow(abs(a'), []);
end


%   ============================================================================
%   The R=2 problem
%   ============================================================================

R       =   [1,2];
kernel  =   [3, 2];

mask    =   false(32,96,96);
mask(:,:,1:2:end)   =   true;

data    =   truth.*mask;

recon   =   grappa(data, truth, R, kernel);
img_combined  = zeros(size(recon));

show_quad(data, recon, 'R=2');
snr_r2 = snr(truth, truth - recon);
mse_r2 = mse(truth, recon);

snr_arr = zeros(32, 1);
mse_arr = zeros(32, 1);
    figure;
for i = 1:32
    temp = truth(i, :, :);

    subplot(1, 2, 1);
    imshow(squeeze(abs(ifft2(ifftshift(temp)))), []);
    subplot(1, 2, 2);
    imshow(squeeze(abs(img_combined)), []);
    snr_arr(i) = snr(abs(temp), temp - img_combined);
    mse_arr(i) = mse(abs(temp), img_combined);
end

disp(snr_arr);
disp(mse_arr);
%% 



%   ============================================================================
%   The R=3 problem
%   ============================================================================

R       =   [1,3];
kernel  =   [3,4];

mask    =   false(32,96,96);
mask(:,:,1:3:end)   =   true;

data    =   truth.*mask;

recon   =   grappa(data, calib, R, kernel);

show_quad(data, recon, 'R=3');


%   ============================================================================
%   The R=6 problem
%   ============================================================================

R       =   [1,6];
kernel  =   [3,2];

mask    =   false(32,96,96);
mask(:,:,1:6:end)   =   true;

data    =   truth.*mask;

recon   =   grappa(data, calib, R, kernel);

show_quad(data, recon, 'R=6');


%   ============================================================================
%   The noisy R=6 problem
%   ============================================================================

R       =   [1,6];
kernel  =   [3,2];

mask    =   false(32,96,96);
mask(:,:,1:6:end)   =   true;

noise   =   1E-6*(randn(size(mask)) + 1j*randn(size(mask)));
data    =   (truth + noise).*mask;

recon   =   grappa(data, calib, R, kernel);

show_quad(data, recon, 'R=6 with noise');