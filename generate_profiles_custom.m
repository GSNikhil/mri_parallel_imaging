function [profiles] = generate_profiles_custom(n_profiles, Nx, Ny, Nz, type, path, display)

profiles = zeros(Nx, Ny, Nz, n_profiles);
img_profiles = zeros(Nx, Ny, Nz, n_profiles); % Only for display

rng(442); % Random Seed 
span_x = linspace(-Nx/2, Nx/2 - 1, Nx);
span_y = linspace(-Ny/2, Ny/2 - 1, Ny);
[x, y] = meshgrid(span_x, span_y);
p = phantom('Modified Shepp-Logan', Nx);

if type == 0
    weights = randn(n_profiles, 2);
    for i = 1:n_profiles
        profiles(:, :, :, i) = 1 ./ (exp(-(weights(i, 1) * x + weights(i, 2) * y)) + 1);
        img_profiles(:, :, :, i) = p .* profiles(:, :, :, i);
    end
elseif type == 1
    centers = randi([-Nx/2, Nx/2 - 1], n_profiles, 2);
    %var = randi([0, 20], n_profiles, 1);
    var = zeros(16, 1);
    for i = 1:n_profiles
        var(i) = 500;
        profile = exp(-((x - centers(i, 1)).^2 + (y - centers(i, 2)).^2) / (2 * var(i))) / (2 * pi * var(i));
        profile = (profile - min(profile(:))) / (max(profile(:)) - min(profile(:)));
        img_profiles(:, :, :, i) = p .* profile;
        profiles(:, :, :, i) = profile;
    end
elseif type == 2
    
end

% Do FFT and save them for SENSE
ffts = zeros(Nx, Ny, Nz, n_profiles);
for i = 1:n_profiles
    ffts(:, :, :, i) = (fft2(img_profiles(:, :, :, i)));
end

% plot all profiles on phantom
if display
    f=figure;
    for i = 1:n_profiles
        subplot(n_profiles / 4, 4, i);
        imshow(img_profiles(:, :, :, i));
    end
    sgtitle("Images with Profiles");
    file_name = "plots/image_profiles_" + n_profiles + "_" +  type + ".png";
    print(f,'-dpng','-r0',file_name);
    
    % plot all profiles
    f=figure;
    for i = 1:n_profiles
        subplot(n_profiles / 4, 4, i);
        imshow(profiles(:, :, :, i));
    end
    sgtitle("Profiles");
    file_name = "plots/profiles_" + n_profiles + "_" +  type + ".png";
    print(f,'-dpng','-r0',file_name);

    f=figure;
    for i = 1:n_profiles
        subplot(n_profiles / 4, 4, i);
        imshow(log10(abs(ffts(:, :, :, i))), []);
    end
    sgtitle("FFTs");
    file_name = "plots/ffts_" + n_profiles + "_" +  type + ".png";
    print(f,'-dpng','-r0',file_name);

end

save(path, "profiles", "img_profiles", "ffts");
end