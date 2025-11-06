%folder = "/home.local/roetzer/C_N/FinalConfig_Nomove_mucilageCN100_PaperReady_spreadNoDivide_leak/";
folder = "/home.local/roetzer/C_N/FinalConfig_paperReady_referenceSetting_exudationDays/" % 

parameters.K_Cliquid =  3.9232* 10^-4;
parameters.minConC_B = 0.0132/8; 
parameters.initConC_B = 0.0539;
parameters.maxConcC_B = 0.3168;

output_file_vtk = folder + "paperBilder/C_BVector/";
mkdir(output_file_vtk)  

output_file_vtk = folder + "paperBilder/C_N_SVector/";
mkdir(output_file_vtk) 

output_file_vtk = folder + "paperBilder/C_SVector/";
mkdir(output_file_vtk) 

output_file_vtk = folder + "paperBilder/UVector/";
mkdir(output_file_vtk) 

output_file_vtk = folder + "paperBilder/N_SVector/";
mkdir(output_file_vtk) 

output_file_vtk = folder + "paperBilder/solution/";
mkdir(output_file_vtk) 

output_file_vtk = folder + "paperBilder/edgeChargeVector/";
mkdir(output_file_vtk)  

output_file_vtk = folder + "paperBilder/C_NVector/";
mkdir(output_file_vtk) 

output_file_vtk = folder + "paperBilder/test/";
mkdir(output_file_vtk) 

for i = 455:5:500
    data = load(folder +"config." + string(i) + ".mat");
    value2vector = ones(data.g.numT, 1);
   % carbon = data.C_BVector + data.C_SVector + C_MNVector + C_PMNVector;
   % nitrogen = data.C_BVector/10 + data.C_SVector + C_MNVector/10 + C_PMNVector/100;
    
    output_file_vtk = folder + "paperBilder/C_BVector/";
   visualizeDataSub(data.g, data.C_BVector, 'C_BVector', 'C_BVector', i,char(output_file_vtk));

     output_file_vtk = folder + "paperBilder/C_N_SVector/";
     visualizeDataSub(data.g, data.C_SVector./data.N_SVector, 'C_N_SVector', 'C_N_SVector', i,char(output_file_vtk));
    % 
    % output_file_vtk = folder + "paperBilder/C_SVector/";
    % visualizeDataSub(data.g, data.C_SVector, 'C_SVector', 'C_SVector', i,char(output_file_vtk));
    % 
     output_file_vtk = folder + "paperBilder/UVector/";
     UVector =((data.C_SVector)./(data.C_SVector + value2vector .* parameters.K_Cliquid));
     %UVector =data.C_SVector;
   visualizeDataSub(data.g, UVector, 'UVector', 'UVector', i,char(output_file_vtk));
    % 
    % output_file_vtk = folder + "paperBilder/N_SVector/";
    % visualizeDataSub(data.g, data.N_SVector, 'N_SVector', 'N_SVector', i,char(output_file_vtk));
    % 
    % output_file_vtk = folder + "paperBilder/C_N_SVector/";
    % visualizeDataSub(data.g, data.C_SVector./data.N_SVector, 'C_N_SVector', 'C_N_SVector', i,char(output_file_vtk));
    
     output_file_vtk = folder + "paperBilder/solution/";
     MB_Vector_vis = data.MB_Vector > ((parameters.maxConcC_B - parameters.minConC_B)/2);
     MN_Vector_vis = (data.MB_Vector > 0 & data.MB_Vector <= (parameters.maxConcC_B - parameters.minConC_B)/2);
     visualizeDataSub(data.g, data.bulkVector + data.POMVector + data.MNVector + MN_Vector_vis *3 + MB_Vector_vis *4  , 'cellType', 'solu', i,char(output_file_vtk));
    % 
    %  output_file_vtk = folder + "paperBilder/edgeChargeVector/";
    % visualizeDataEdges(data.g, data.edgeChargeVector, 'memoryEdges', 'edgeChargeVector', i, 2, char(output_file_vtk));


    %output_file_vtk = folder + "paperBilder/C_NVector/";
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

     % output_file_vtk = folder + "paperBilder/test/";
     % C_N = data.C_SVector./data.N_SVector;
     % C_N(:) = 0;
     % C_N(1) = 100;
     % visualizeDataSub(data.g, C_N, 'C_N', 'C_N', i,char(output_file_vtk));
     

end