%   examples.m
%   mchiew@fmrib.ox.ac.uk

%   Load example data
R_range = [1, 2, 3, 4, 6];
coil_type_range = [0, 1];
coil_count_range = [4, 8, 16];
display_profiles = false;

Nx = 96;
Ny = 96;
Nz = 1;

p = phantom('Modified Shepp-Logan', Nx);

grappa_imgs = {};
sense_imgs = {};

for type = coil_type_range
    for n_coil = coil_count_range
        path = strcat('data/profile' , num2str(type) , '.mat');
        generate_profiles_custom(n_coil, Nx, Ny, Nz, type, path, display_profiles);
    
        input   =   matfile(path);
        truth   =   squeeze(input.ffts);
        calib   =   squeeze(input.ffts);
        
        truth = permute(truth, [3 2 1]);
        calib = permute(calib, [3 2 1]);
    
        if display_profiles
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
        end
    
        for skip_factor = R_range
            % Grappa Begin
            R       =   [1, skip_factor];
            kernel  =   [3, 2];
    
            mask    =   false(n_coil,Nx,Ny);
            mask(:,:,1:skip_factor:end)   =   true;
    
            data    =   truth.*mask;
    
            recon   =   grappa(data, truth, R, kernel);
            plot_title = ['R=', num2str(skip_factor), ' #Coils=', num2str(n_coil), ' Type=', num2str(type)];
            %show_quad(data, recon, plot_title);
   
            img_combined  = zeros(size(recon));
            for i = 1:n_coil
                img_combined(i, :, :) = (ifft2(squeeze(recon(i, :, :)))); 
            end
    
            img_combined = sum(img_combined.^2, 1).^0.5;
            img_combined = abs(squeeze(img_combined))';
            
            [mse_grappa, snr_grappa] = compute_metrics(p, img_combined);
            

            %plot_title = [plot_title, "Grappa", "MSE=", mse_grappa, "SNR=", snr_grappa];
            
            % Grappa End

            % Sense Begin
            sense_recon = sense(input.ffts, skip_factor);
            [mse_sense, snr_sense] = compute_metrics(p, abs(sense_recon));
            
            %plot_title = [plot_title, ", Sense", "MSE=", mse_sense, "SNR=", snr_sense];
            % Sense End
            
            %disp(plot_title);

            grappa_imgs = [grappa_imgs, img_combined];
            sense_imgs = [sense_imgs, abs(sense_recon)];
            
            f=figure;
            subplot(1, 3, 1);
            imshow(p, []);
            title("Original Image");

            subplot(1, 3, 2);         
            imshow(img_combined, []);
            grappa_title = ["Grappa Recon", "MSE: " + num2str(mse_grappa), "SNR: " + num2str(snr_grappa)]; 
            title(grappa_title);

            subplot(1, 3, 3);
            imshow(normalize_max_min(abs(sense_recon)), []);
            sense_title = ["Sense Recon" , "MSE: " + num2str(mse_sense), "SNR: " + num2str(snr_sense)]; 
            title(sense_title);

            sgtitle(plot_title);
            file_name = "plots/" + n_coil + "_" +  type + "_" + skip_factor + ".png";
            %saveas(f, file_name);
            print(f,'-dpng','-r0',file_name);
            
            % Type - nCoils - R - Grappa MSE - Grappa SNR - Sense MSE -
            % Sense SNR
            stat_arr = [type, n_coil, skip_factor, mse_grappa, snr_grappa, mse_sense, snr_sense];
            disp(stat_arr);
        end
    end
end
plot_comparisions(0, grappa_imgs, "Grappa Type 0", 0);
plot_comparisions(1, grappa_imgs, "Grappa Type 1", 0);
plot_comparisions(0, sense_imgs, "Sense Type 0", 1);
plot_comparisions(1, sense_imgs, "Sense Type 1", 1);

function plot_comparisions(type, imgs, plot_title, normalize)
    f=figure;
    idx = 15 * type + 1;
    for n_coil = [4, 8, 16]
        for skip_factor = [1, 2, 3, 4, 6]
            subplot(3, 5, idx - 15 * type);
            img = imgs{idx};
            if normalize
                img = normalize_max_min(img);
            end
            imshow(img, []);
            title("R = " + num2str(skip_factor) + " nCoils = " + num2str(n_coil));
            idx = idx + 1;
        end
    end
    sgtitle(plot_title);
    file_name = "plots/" + plot_title + "_all.png";
    print(f,'-dpng','-r0',file_name);
end

function [mse, snr] = compute_metrics(true_img, recon_img)
    %recon_img = normalize_max_min(recon_img);
    
    %{
    figure;
    subplot(1, 3, 1);
    imshow(true_img, []);
    subplot(1, 3, 2);
    imshow(recon_img, []);
    subplot(1, 3, 3);
    imshow(recon_img - true_img, []);
    %}
    
    error = recon_img - true_img;
    error_std = var(error(:));
    
    squared_error = (error .^ 2);
    mse = mean(squared_error(:));
    snr = mean(mean(recon_img)) / error_std;
end

function [img_n] = normalize_max_min(img)
    maxi = max(img);
    mini = min(img);
    img_n = (img - mini) ./ (maxi - mini);
end