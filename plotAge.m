function plotAge()
close all
root = "/home.local/roetzer/C_N/FinalConfig_O2_2/";
%root = "FinalConfig/"  
config_file = 'config.50.mat';


file =  root + config_file
data = load(file);
image = reshape(data.POMageVector, [250, 250]);  
image(image > 0) =  ceil(image(image > 0) /24) ;
image = image .* reshape(data.MNVector, [250, 250]);
bulkImage = reshape(data.bulkVector, [250, 250]);
%figure
%imshow(image, [])
figure
image2 = image;
%image2(image == 0) = bulkImage(image == 0) * 100;
maxi = ceil(max(image2,[],'all'))
grey_value = maxi +2
inverse = ones([250,250]) .* maxi;
%image2 = inverse - image2;
%image2(image2 == maxi) = 0;
%imshow(image2, [])
%colormap summer    
mineralPhase = bulkImage;
mineralPhase(image > 0) = 0;
mineralPhase = grey_value *  mineralPhase;


days_matrix = image +  mineralPhase ;
maxi = ceil(max(image2,[],'all'))
imagesc(days_matrix);

% Add color bar
colorbar_handle = colorbar;

% Label the colorbar
ylabel(colorbar_handle, 'Days');
% Optionally, add labels and title
%xlabel('X-axis label');
%ylabel('Y-axis label');
title('Age of Nercomass - constant Decay');

% Adjust color map if necessary (optional)
%colormap(jet); % You can choose other colormaps like 'parula', 'hot', 'cool', etc.

% Define a custom colormap going from dark green to light green
% and make zero values white
custom_colormap = [
    0.0 0.2 0.0;  % dark green
    0.2 0.4 0.2;  % darker green
    0.4 0.6 0.4;  % medium green
    0.6 0.8 0.6;  % lighter green
    0.8 1.0 0.8;  % light green
];

% Interpolate the colormap to get a smooth gradient
num_colors = 256; % Number of colors in the colormap
custom_colormap = interp1(linspace(0, 1, size(custom_colormap, 1)), custom_colormap, linspace(0, 1, num_colors));

% Insert white for zero values at the beginning of the colormap
custom_colormap = [[1 1 1]; custom_colormap];

% % Determine the normalized index for the grey value in the colormap
% min_val = min(days_matrix(:));
% max_val = max(days_matrix(:));
% grey_index = round((grey_value - min_val) / (max_val - min_val) * (num_colors - 1)) + 1;
% 
% % Ensure the grey value index is within the valid range
% grey_index = min(max(grey_index, 1), num_colors + 1); % +1 because we added white
% 
% % Set the grey color at the determined index
% custom_colormap(grey_index, :) = [0.5 0.5 0.5]; % Grey color
custom_colormap = [custom_colormap; [0.5 0.5 0.5]];
% Apply the custom colormap
colormap(custom_colormap);


% Adjust the color limits to include all data
caxis([min(days_matrix(:)), max(days_matrix(:))]);
end