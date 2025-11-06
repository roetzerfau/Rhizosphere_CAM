%folder_input = "/home.other/fauamlx5g/roetzer/C_N/FinalConfig_paperReady_referenceSetting_exudationDays/"; % 
folder_input = "/home.other/fauamlx5g/roetzer/C_N/FinalConfig_paperReady_referenceSetting_CN10/"; % 
folder_input = "/home.local/roetzer/C_N/FinalConfig_paperSimulation_CN40_referenceSetting/"; % 
folder_input = "/home.local/roetzer/C_N/FinalConfig_paperSimulation_referenceSetting_CN40_without_particle_Movement/"; % 
%paperReady_referenceSetting_CN40_without_particle_Movement
folder_output = pwd+"/finalpics/paperSimulation_CN40_referenceSetting";
folder_output = pwd+"/finalpics/paperSimulation_CN40_referenceSetting_without_particle_Movement";
parameters.K_Cliquid =  3.9232* 10^-4;
parameters.minConC_B = 0.0132/8; 
parameters.initConC_B = 0.0539;
parameters.maxConcC_B = 0.3168;

output_file_vtk = folder_output + "paperBilder/C_BVector/";
%mkdir(output_file_vtk)  

output_file_vtk = folder_output + "paperBilder/C_N_SVector/";
mkdir(output_file_vtk) 

output_file_vtk = folder_output + "paperBilder/C_SVector/";
%mkdir(output_file_vtk) 

output_file_vtk = folder_output + "paperBilder/UVector/";
%mkdir(output_file_vtk) 

output_file_vtk = folder_output + "paperBilder/N_SVector/";
%mkdir(output_file_vtk) 

output_file_vtk = folder_output + "paperBilder/solution/";
%mkdir(output_file_vtk) 

output_file_vtk = folder_output + "paperBilder/edgeChargeVector/";
%mkdir(output_file_vtk)  

output_file_vtk = folder_output + "paperBilder/C_NVector/";
%mkdir(output_file_vtk) 

output_file_vtk = folder_output + "paperBilder/test/";
%mkdir(output_file_vtk) 

colors = struct(... % used from here ->
    'pom', hex2rgb('#84CC83'), ...
    'new_microbial_necromass', hex2rgb('#CC9922'), ...
    'solid', hex2rgb('#808080'), ...
    'c_n_s', hex2rgb('#4287f5'), ...
    'dissolved_substrate', hex2rgb('#E64640'), ... 
    'microbial_biomass', hex2rgb('#AAFF00'), ...
    'old_microbial_necromass', hex2rgb('#AA5500') ...
);
colors_gradient = struct(...
    'u_start', [0 0 1],... % hex2rgb('#395DBF'),... [0 0 1]
    'u_end',  [1 1 1],... % hex2rgb('#A8323E')...
    'c_n_start', [1 0 0],... % hex2rgb('#395DBF'),... [0 0 1]
    'c_n_end',  [1 1 1]... % hex2rgb('#A8323E')...
);

isVisible = 'on';

for i = [125:1:135]%455:5:455%500  195:5:230
    data = load(folder_input +"config." + string(i) + ".mat");
    value2vector = ones(data.g.numT, 1);
    % carbon = data.C_BVector + data.C_SVector + C_MNVector + C_PMNVector;
    % nitrogen = data.C_BVector/10 + data.C_SVector + C_MNVector/10 + C_PMNVector/100;
    
    output_file_vtk = folder_output + "paperBilder/C_BVector/";
    %visualizeDataSub(data.g, data.C_BVector, 'C_BVector', 'C_BVector', i,char(output_file_vtk));
    
    
    output_file_vtk = folder_output + "paperBilder/C_N_SVector/";
    %!!! visualizeDataSub(data.g, data.C_SVector./data.N_SVector, 'C_N_SVector', 'C_N_SVector', i,char(output_file_vtk));
    % 
    % output_file_vtk = folder_output + "paperBilder/C_SVector/";
    % visualizeDataSub(data.g, data.C_SVector, 'C_SVector', 'C_SVector', i,char(output_file_vtk));
    % 
    %!!! output_file_vtk = folder_output + "paperBilder/UVector/";
    UVector =((data.C_SVector)./(data.C_SVector + value2vector .* parameters.K_Cliquid));
    %UVector =data.C_SVector;
    % visualizeDataSub(data.g, UVector, 'UVector', 'UVector', i,char(output_file_vtk));
    % 
    % output_file_vtk = folder_output + "paperBilder/N_SVector/";
    % visualizeDataSub(data.g, data.N_SVector, 'N_SVector', 'N_SVector', i,char(output_file_vtk));
    % 
    % output_file_vtk = folder_output + "paperBilder/C_N_SVector/";
    % visualizeDataSub(data.g, data.C_SVector./data.N_SVector, 'C_N_SVector', 'C_N_SVector', i,char(output_file_vtk));
    
    %!!! output_file_vtk = folder_output + "paperBilder/solution/";
    MB_Vector_vis = data.MB_Vector > ((parameters.maxConcC_B - parameters.minConC_B)/2);
    MN_Vector_vis = (data.MB_Vector > 0 & data.MB_Vector <= (parameters.maxConcC_B - parameters.minConC_B)/2);
    %!!!! visualizeDataSub(data.g, data.bulkVector + data.POMVector + data.MNVector + MN_Vector_vis *3 + MB_Vector_vis *4  , 'cellType', 'solu', i,char(output_file_vtk));
    %imshow(reshape(data.POMVector,[data.g.NX,data.g.NX]),[])
    mn_new = zeros(size(data.MNVector));
    mn_new(data.POMageVector~=max(data.POMageVector)) = data.MNVector(data.POMageVector~=max(data.POMageVector));
    mn_old = zeros(size(data.MNVector));
    mn_old(data.POMageVector==max(data.POMageVector)) = data.MNVector(data.POMageVector==max(data.POMageVector));
    %POMageVector(index(MNVector)) > 1 -> new
    
    layers = { %for values 1 or 0, one color
        reshape(data.POMVector,[data.g.NX,data.g.NX]), colors.pom;
        reshape(mn_new,[data.g.NX,data.g.NX]), colors.new_microbial_necromass;
        reshape(mn_old,[data.g.NX,data.g.NX]), colors.old_microbial_necromass;
        reshape(data.bulkVector - data.POMVector - data.MNVector,[data.g.NX,data.g.NX]), colors.solid;
        reshape(data.C_SVector,[data.g.NX,data.g.NX]), colors.dissolved_substrate;
    };

    % Color Grading 1 Color for values >0, adjustable with conc_mult:
    % double between [0,1] 0: highest concentration is black, 1:all concentrations equal
    conc_mult=1; 
    layers_conc = { %for values as gradients one color
        reshape(data.C_BVector,[data.g.NX,data.g.NX]), colors.microbial_biomass;
    };

    color_var = '';
    %color_var = '2colors_';
    %color_var = '1color_';
    output_file_vtk = folder_output + "paperBilder/C_N_SVector/";
    %output_file_vtk = folder_output + "paperBilder/UVector/";
    
    % Color Grading 2 Colors for values between 0,1 equally
    ls_grad_2col_v1 = {
       %reshape(UVector,[data.g.NX,data.g.NX]), colors_gradient.u_start, colors_gradient.u_end;
    };
    
    % Color Grading 2 Colors for values >0, mid point and boundary adjustable
    ls_grad_2col_v2 = {%data.C_SVector./data.N_SVector
        %reshape(UVector,[data.g.NX,data.g.NX]), colors_gradient.u_start, colors_gradient.u_end;
    };

    %color grading with ounly 1 color
    ls_grad_1col = { %for values as gradients one color
       reshape(data.C_SVector./data.N_SVector,[data.g.NX,data.g.NX]), colors_gradient.c_n_start, 150;
       %reshape(UVector,[data.g.NX,data.g.NX]), colors_gradient.u_start, 1;
    };
    
    rgb_img = ones(data.g.NX, data.g.NX, 3);
    
    
    for j = 1:size(ls_grad_2col_v1,1)
        
        
        value = ls_grad_2col_v1{j,1};
        
        start_color =ls_grad_2col_v1{j,2};
        end_color =ls_grad_2col_v1{j,3};
        max_val = 1;
        min_val = 0;
        
        mask_lower =  (value <max_val/2) & (value >min_val);
        mask_upper = (value >= max_val/2);
        for c = 1:3
            temp = rgb_img(:,:,c);
            lowest_color_fac =0.1;
            
            val_lower = (1 - 2*value(mask_lower)/(max_val));
            val_lower_withMin = lowest_color_fac + val_lower*(1-lowest_color_fac);
            temp(mask_lower) = 1 - val_lower_withMin * (1-start_color(c));
            
            val_upper = 2*(value(mask_upper)/(max_val)-0.5)+1;
            val_upper_withMin = lowest_color_fac + val_upper*(1-lowest_color_fac);
            temp(mask_upper) = 1- val_upper_withMin * (1-end_color(c));
            
            rgb_img(:,:,c) = temp;
        end
    end
    
    for j = 1:size(ls_grad_2col_v2,1)
        
        
        value = ls_grad_2col_v2{j,1};
        
        start_color =ls_grad_2col_v2{j,2};
        end_color =ls_grad_2col_v2{j,3};
        max_val = 100;
        min_val = 0;
        
        mask_lower =  (value <max_val/2) & (value >min_val);
        mask_upper = (value >= max_val/2) & (value < max_val);
        mask_upper_bound = (value >= max_val);
        for c = 1:3
            temp = rgb_img(:,:,c);
            lowest_color_fac =0.1;
            
            val_lower = (1 - 2*value(mask_lower)/(2*max_val));
            val_lower_withMin = lowest_color_fac + val_lower*(1-lowest_color_fac);
            temp(mask_lower) = 1 - val_lower_withMin * (1-start_color(c));
            
            val_upper = 2*(value(mask_upper)/(2*max_val)-0.5)+1;
            val_upper_withMin = lowest_color_fac + val_upper*(1-lowest_color_fac);
            temp(mask_upper) = 1- val_upper_withMin * (1-end_color(c));
            
            temp(mask_upper_bound) = end_color(c);
            
            rgb_img(:,:,c) = temp;
            
        end
    end

    for j = 1:size(ls_grad_1col,1)
        value = ls_grad_1col{j,1};
        color = ls_grad_1col{j,2};
        upper_bound = ls_grad_1col{j,3};

        mask =  value <= upper_bound;
        mask_bound = value > upper_bound;
        max_el = max(value(:));
        for c = 1:3
            lowest_color_fac =0.1;
            
            temp = rgb_img(:,:,c);
            val = (value(mask)/upper_bound);
            val_lower_withMin = lowest_color_fac + val*(1-lowest_color_fac);
            temp(mask) =  1-val_lower_withMin * (1-(color(c)));
            temp(mask_bound) = color(c);
            rgb_img(:,:,c)=temp;
        end
    end

    for j = 1:size(layers,1)
        mask = layers{j,1} == 1;
        for c = 1:3
            temp = rgb_img(:,:,c);
            temp(mask) = layers{j,2}(c);
            rgb_img(:,:,c)=temp;
        end
    end

    
    for j = 1:size(layers_conc,1)
        mask = layers_conc{j,1} ~= 0;
        max_el = max(layers_conc{j,1}(:));
        for c = 1:3
            temp = rgb_img(:,:,c);
            temp(mask) = (conc_mult+(1-layers_conc{j,1}(mask)/max_el)*(1-conc_mult)) * layers_conc{j,2}(c);
            rgb_img(:,:,c)=temp;
        end
    end


    
    
    
    figure1 = figure('visible', isVisible);
    imshow(rgb_img)
    %imwrite(rgb_img, fullfile(output_file_vtk, "solu2" + ".png"))
    
    set(figure1, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    pause(1)
    exportgraphics(figure1, fullfile(output_file_vtk, color_var +"day" +string(i)+ ".png"),'Resolution',400);
    
    %savefig(fullfile(output_file_vtk, "solu"))
    % 
    %  output_file_vtk = folder_output + "paperBilder/edgeChargeVector/";
    % visualizeDataEdges(data.g, data.edgeChargeVector, 'memoryEdges', 'edgeChargeVector', i, 2, char(output_file_vtk));
    
    %solution = soikparticles + pom + necro + biomass + dissovled 
    %bild 
    %output_file_vtk = folder_output + "paperBilder/C_NVector/";
    %C_N = data.C_SVector./data.N_SVector;
    %%C_N(isnan(C_N)) = 0;
    %%C_N = C_N + ((data.POMVector -data.MNVector) > 0) * 100 + (data.MB_Vector > 0) * 10 + (data.MNVector > 0) * 10; 
    %%C_N((data.bulkVector - data.POMVector) > 0) = NaN;
    
    %%imshow(reshape(C_N,[data.g.NX,data.g.NX]),[])
    %visualizeDataSub(data.g, C_N, 'C/N', 'CN', i,char(output_file_vtk));
    
    % % Define color map for the range 0-100
    % colormapRange = jet(256); % Use a colormap like 'jet'
    % 
    % % Extend the colormap for out-of-range values
    % outOfRangeColor = [200, 200, 200] / 255; % White color for out-of-range
    % colormapFull = [colormapRange; outOfRangeColor];
    % 
    % colorIndices = zeros(size(C_N));
    % colorIndices(C_N >= 0 & C_N <= 100) = round((C_N(C_N>= 0 & C_N <= 100) / 100) * 255) + 1; % Map to colormap
    % colorIndices(isnan(C_N)) = 257; % Last color (out-of-range high)
    % 
    % figure;
    % imagesc(reshape(colorIndices,[data.g.NX,data.g.NX]));
    % colormap(colormapFull);
    % %colorbar('Ticks', [1, 2, 130, 258], ...
    % %         'TickLabels', {'< 0', '0', '50', '> 100'});
    % axis image;
    % axis off;
    
    % output_file_vtk = folder_output + "paperBilder/test/";
    % C_N = data.C_SVector./data.N_SVector;
    % C_N(:) = 0;
    % C_N(1) = 100;
    % visualizeDataSub(data.g, C_N, 'C_N', 'C_N', i,char(output_file_vtk));
         

end
system("mogrify -trim " + output_file_vtk + "*.png")
close all