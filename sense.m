function [recon_img] = sense(raw, R)
    [Nx,Ny,Nz,Nc] = size(raw);
    %show_grid(log(abs(raw)), [], jet);
    img = (ifft2(raw));    
    %show_grid(abs(img),[], gray);

    %%
    %Computing square root of sum of squares of images
    
    img_combined  = sqrt(sum(abs(img).^2,4));
    %show_img(abs(img_combined),[]);
    S_0 = img./img_combined;
    %show_grid(abs(S_0),[],jet);
    %% 
    %Smoothing
    kernel = ones(3)/3^2;
    S_1 = zeros(Nx,Ny,Nz,Nc);
    for i = 1:Nc
        S_1(:,:,1,i) = conv2(S_0(:,:,1,i),kernel,'same');
    end
    %show_grid(abs(S_1),[],jet);
    %Masking sensitivities
    thresh = 0.0001*max(abs(img_combined(:)));
    mask = abs(img_combined) > thresh;
    %%
    S_2 = S_1.*mask;
    %show_grid(abs(S_2),[],jet);
    %% Use sensivities 
    
    img_combined_opt = sum(img.*conj(S_2),4)./sum(S_2.*conj(S_2),4);
    %show_img(abs(img_combined_opt),[],gray);
    
    %%
    raw_R3 = raw;
    for i = 1:R-1
        raw_R3(i:R:end,:,:,:) = 0;
    end
    img_R3 = (ifft2(raw_R3));
    %show_grid(abs(img_R3),[],gray);
    recon_img = sense_helper(img_R3,S_2,R);
    %show_img(abs(recon_img),[],gray);
end

function out = sense_helper(input,sens,R)
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