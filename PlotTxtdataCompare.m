function PlotTxtdataCompare(files, description, names_appendix)
close all
imageFolder = "bilder/CN/paperReady/";
%imageFolder = "final_pics/";
outputfolder = imageFolder + description + "/";
if ~exist(outputfolder, 'dir')  % Check if the folder exists
    mkdir(outputfolder);        % Create the folder
    disp(['Folder created at: ', outputfolder]);
else
    disp(['Folder already exists at: ', outputfolder]);
end
%ANNZ / SNNZ    * (sum(C_A)/ANNZ)*  (C_SNNZ/sum(C_Solid)  
%flÃ¤che anteil/ flÃ¤che soil   *  anteil densitiy/solid density 
porosity = 0.411;
soil_particleNNZ = 250 * 250 * (1-porosity);
soilParticleDensity = 2.36;
isNNZ = false;
startValue = 5;
soilfactor = ( soil_particleNNZ * soilParticleDensity/1000 );

titel = description;
percentage_plot = 50;%;%25;
isVisible = 'on';
average_factor = 1;
if files.contains("2years")
    last_value = 2000;
else
last_value = Inf;
end
% DOC ist zu hoch

ylim_C_stacked = [0 0.2];
ylim_N_stacked = [0 4*10^-3];

ylim_C = [10^-8 10^-0];
ylim_N = [10^-8 10^0 ];



%colorCode = get(gca,'colororder');
colors = struct(...
    'microbial_biomass', hex2rgb('#AAFF00'), ...
    'new_microbial_necromass', hex2rgb('#CC9922'), ...
    'old_microbial_necromass', hex2rgb('#AA5500'), ...
    'sum_microbial_necromass', hex2rgb('#FF8800'), ...
    'pom', hex2rgb('#84CC83'), ...
    'pom_added', hex2rgb('#B4E0B4'), ...
    'dissolved_substrate', hex2rgb('#E64640'), ...
    'leakage_co2', hex2rgb('#000000'), ...
    'leakage_co2_text', hex2rgb('#FFFFFF'), ... % for barplot text
    'r_o', hex2rgb('#79C7F7'), ...
    'r', hex2rgb('#000000'), ...  % also r_leak
    'cue', hex2rgb('#000000'), ...
    'c_n', hex2rgb('#000000'), ...
    'unknown', hex2rgb('#FC08EC') ...  
);

fontsizes = struct(...
    'xlabel', 16, ...
    'ylabel', 16, ... 
    'xaxis', 16, ...
    'yaxis', 16, ... 
    'legend', 18, ...
    'title', 18, ...
    'bartext', 12, ...
    'bartext_arrow', 12 ...
);
title_on = false;
tick_dir = 'out'; %or 'in'
    
close all
%% Plots 
linewidth = 2;
markers = ["-", "--", ":","-."];
if numel(files) <= numel(markers)
%% figure 5
figure5 = figure('visible', isVisible) 
names = {};
for file = 1:numel(files)
appendix = files(file);
root = "/home.local/roetzer/C_N/txtdata_" + appendix +"/"
%root = "/home.other/fauamlx5g/roetzer/C_N/txtdata_" + appendix +"/"
    
%% carbon data
fileID = fopen(root + 'C_BVector.txt','r');
formatSpec = '%f %f';
C_B = fscanf(fileID,formatSpec,  [2 last_value])';
y_C_B = C_B(startValue:end,2)/average_factor;
x_C_B = C_B(startValue:end,1);

if(isNNZ)
fileID = fopen(root + 'N_BVectorNNZ.txt','r');
formatSpec = '%f %f';
C_B = fscanf(fileID,formatSpec,  [2 last_value])';
y_C_BNNZ = C_B(startValue:end,2);
x_C_BNNZ = C_B(startValue:end,1);
end

fileID = fopen(root + 'C_SVector.txt','r');
formatSpec = '%f %f';
C_S = fscanf(fileID,formatSpec,  [2 last_value])';
y_C_S = C_S(startValue:end,2)/average_factor;
x_C_S = C_S(startValue:end,1);

if(isNNZ)
fileID = fopen(root + 'C_SVectorNNZ.txt','r');
formatSpec = '%f %f';
C_S = fscanf(fileID,formatSpec,  [2 last_value])';
y_C_SNNZ = C_S(startValue:end,2);
x_C_SNNZ = C_S(startValue:end,1);
end

fileID = fopen( root + 'C_MNVector.txt','r');
formatSpec = '%f %f';
C_MN = fscanf(fileID,formatSpec,  [2 last_value])';
y_C_MN = C_MN(startValue:end,2)/average_factor;
x_C_MN = C_MN(startValue:end,1);


if(isNNZ)
fileID = fopen( root + 'C_MNVectorNNZ.txt','r');
formatSpec = '%f %f';
C_MN = fscanf(fileID,formatSpec,  [2 last_value])';
y_C_MNNNZ = C_MN(startValue:end,2);
x_C_MNNNZ = C_MN(startValue:end,1);
end

fileID = fopen( root + 'C_MNVector_old.txt','r');
formatSpec = '%f %f';
C_MN_old = fscanf(fileID,formatSpec,  [2 last_value])';
y_C_MN_old = C_MN_old(startValue:end,2)/average_factor;
x_C_MN_old = C_MN_old(startValue:end,1);

if(isNNZ)
fileID = fopen( root + 'C_MNVector_oldNNZ.txt','r');
formatSpec = '%f %f';
C_MN_old = fscanf(fileID,formatSpec,  [2 last_value])';
y_C_MN_oldNNZ = C_MN_old(startValue:end,2);
x_C_MN_oldNNZ = C_MN_old(startValue:end,1);
end

fileID = fopen( root + 'C_MNVector_new.txt','r');
formatSpec = '%f %f';
C_MN_new = fscanf(fileID,formatSpec,  [2 last_value])';
y_C_MN_new = C_MN_new(startValue:end,2)/average_factor;
x_C_MN_new = C_MN_new(startValue:end,1);

if(isNNZ)
fileID = fopen( root + 'C_MNVector_newNNZ.txt','r');
formatSpec = '%f %f';
C_MN_new = fscanf(fileID,formatSpec,  [2 last_value])';
y_C_MN_newNNZ = C_MN_new(startValue:end,2);
x_C_MN_newNNZ = C_MN_new(startValue:end,1);
end

fileID = fopen( root + 'C_POMconcVector.txt','r');
formatSpec = '%f %f';
C_POM = fscanf(fileID,formatSpec,  [2 last_value])';
y_C_POM = C_POM(startValue:end,2)/average_factor;
x_C_POM = C_POM(startValue:end,1);

if(isNNZ)
fileID = fopen( root + 'C_POMconcVectorNNZ.txt','r');
formatSpec = '%f %f';
C_POM = fscanf(fileID,formatSpec, [2 last_value])';
y_C_POMNNZ = C_POM(startValue:end,2);
x_C_POMNNZ = C_POM(startValue:end,1);
end

fileID = fopen( root + 'CO2Vector.txt','r');
formatSpec = '%f %f';
CO2 = fscanf(fileID,formatSpec,  [2 last_value])';
y_CO2 = CO2(startValue:end,2)/average_factor;
x_CO2 = CO2(startValue:end,1);

if(isNNZ)
fileID = fopen( root + 'CO2VectorNNZ.txt','r');
formatSpec = '%f %f';
CO2 = fscanf(fileID,formatSpec,  [2 last_value])';
y_CO2NNZ = CO2(startValue:end,2);
x_CO2NNZ = CO2(startValue:end,1);
end

fileID = fopen( root + 'CO2Vector_over.txt','r');
formatSpec = '%f %f';
CO2_over = fscanf(fileID,formatSpec,  [2 last_value])';
y_CO2_over = CO2_over(startValue:end,2)/average_factor;
x_CO2_over = CO2_over(startValue:end,1);

if(isNNZ)
fileID = fopen( root + 'CO2Vector_overNNZ.txt','r');
formatSpec = '%f %f';
CO2_over = fscanf(fileID,formatSpec,  [2 last_value])';
y_CO2_overNNZ = CO2_over(startValue:end,2);
x_CO2_overNNZ = CO2_over(startValue:end,1);
end

x = x_C_MN;

POM_added = floor([1:numel(y_C_B)]./10) * 0.083;

plot(x,y_C_B./soilfactor,markers(file),'LineWidth', linewidth, 'Color',colors.microbial_biomass)
hold on
% plot(x,y_C_MN_old./soilfactor,'LineWidth', linewidth)
% hold on
% plot(x,y_C_MN_new./soilfactor,'LineWidth', linewidth)
% hold on
plot(x,y_C_POM./soilfactor,markers(file),'LineWidth', linewidth, 'Color',colors.pom)
hold on
plot(x,y_C_S./soilfactor,markers(file),'LineWidth', linewidth,'Color',colors.dissolved_substrate)
hold on
plot(x,y_CO2./soilfactor,markers(file),'LineWidth', linewidth,'Color',colors.leakage_co2)
hold on
plot(x,(y_C_MN_new + y_C_MN_old)./soilfactor,markers(file),'LineWidth', linewidth,'Color',colors.sum_microbial_necromass)
hold on
qw{file} = plot(x,(y_C_MN_new + y_C_MN_old)./soilfactor * 0,markers(file),'LineWidth', 0.01,'Color',[0 0 0]);
% plot(x,POM_added,'LineWidth', linewidth)
%hold on
%plot(x,y_CO2_over,'LineWidth', linewidth)
name = {'biomass' + names_appendix(file),  ...%'Necromass old'+names_appendix(file), 'Necromass new'+names_appendix(file),...
    'POM'+names_appendix(file), 'dissolved substrate'+names_appendix(file),'CO2'+names_appendix(file), ...
    'necromass total'+names_appendix(file)}; %, 'POM added'+names_appendix(file)};
names = [names,name];



meanC_B = (mean(y_C_B(200:350))- mean(y_C_B(100:199)))./soilfactor


end % files

%legend(names, 'FontSize', 14, 'Location','southeast')
legend([qw{:}],  names_appendix, 'FontSize', fontsizes.legend, 'Location','southeast')
ax = gca
ax.XAxis.FontSize = fontsizes.xaxis
ax.YAxis.FontSize = fontsizes.yaxis
xlabel('days', 'FontSize',fontsizes.xlabel)

ylabel('amount carbon mg C g soil^{-1} (log scale)', 'FontSize', fontsizes.ylabel)
set(gca, 'YScale', 'log');  % Set y-axis to log scale
set(gca,'TickDir',tick_dir);

%set(gca, 'FontSize', 14)

if title_on
    title(titel, 'Interpreter', 'none', 'FontSize',fontsizes.title)
end
set(figure5, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure5, fullfile(outputfolder, "C_log_" +description + ".png"));
savefig(fullfile(outputfolder, "C_log_" + description))
%% figure 6
names = {};
figure6 = figure('visible', isVisible) 
for file = 1:numel(files)
    appendix = files(file);
    %root = "/home.local/roetzer/C_N/txtdata_" +  appendix+"/"
    root = "/home.other/fauamlx5g/roetzer/C_N/txtdata_" +  appendix+"/"
    
%% Nitrogen data
%indices = find(POMageVector > 0 & POMageVector ~= max(POMageVector))
fileID = fopen(root + 'N_BVector.txt','r');
formatSpec = '%f %f';
N_B = fscanf(fileID,formatSpec,  [2 last_value])';
y_N_B = N_B(startValue:end,2)/average_factor;
x_N_B = N_B(startValue:end,1);

if(isNNZ)
fileID = fopen(root + 'N_BVectorNNZ.txt','r');
formatSpec = '%f %f';
N_B = fscanf(fileID,formatSpec,  [2 last_value])';
y_N_BNNZ = N_B(startValue:end,2);
x_N_BNNZ = N_B(startValue:end,1);
end

fileID = fopen(root + 'N_SVector.txt','r');%NNZ
formatSpec = '%f %f';
N_S = fscanf(fileID,formatSpec,  [2 last_value])';
y_N_S = N_S(startValue:end,2)/average_factor;
x_N_S = N_S(startValue:end,1);

if(isNNZ)
fileID = fopen(root + 'N_SVectorNNZ.txt','r');
formatSpec = '%f %f';
N_S = fscanf(fileID,formatSpec,  [2 last_value])';
y_N_SNNZ = N_S(startValue:end,2);
x_N_SNNZ = N_S(startValue:end,1);
end

fileID = fopen( root + 'N_MNVector.txt','r');
formatSpec = '%f %f';
N_MN = fscanf(fileID,formatSpec,  [2 last_value])';
y_N_MN = N_MN(startValue:end,2)/average_factor;
x_N_MN = N_MN(startValue:end,1);

if(isNNZ)
fileID = fopen( root + 'N_MNVectorNNZ.txt','r');
formatSpec = '%f %f';
N_MNNNZ = fscanf(fileID,formatSpec,  [2 last_value])';
y_N_MNNNZ = N_MNNNZ(startValue:end,2);
x_N_MNNNZ = N_MNNNZ(startValue:end,1);
end

fileID = fopen( root + 'N_MNVector_old.txt','r');
formatSpec = '%f %f';
N_MN_old = fscanf(fileID,formatSpec,  [2 last_value])';
y_N_MN_old = N_MN_old(startValue:end,2)/average_factor;
x_N_MN_old = N_MN_old(startValue:end,1);

if(isNNZ)
fileID = fopen( root + 'N_MNVector_oldNNZ.txt','r');
formatSpec = '%f %f';
N_MN_old = fscanf(fileID,formatSpec,  [2 last_value])';
y_N_MN_oldNNZ = N_MN_old(startValue:end,2);
x_N_MN_oldNNZ = N_MN_old(startValue:end,1);
end

fileID = fopen( root + 'N_MNVector_new.txt','r');
formatSpec = '%f %f';
N_MN_new = fscanf(fileID,formatSpec,  [2 last_value])';
y_N_MN_new = N_MN_new(startValue:end,2)/average_factor;
x_N_MN_new = N_MN_new(startValue:end,1);

if(isNNZ)
fileID = fopen( root + 'N_MNVector_newNNZ.txt','r');
formatSpec = '%f %f';
N_MN_new = fscanf(fileID,formatSpec,  [2 last_value])';
y_N_MN_newNNZ = N_MN_new(startValue:end,2);
x_N_MN_newNNZ = N_MN_new(startValue:end,1);
end

fileID = fopen( root + 'N_POMconcVector.txt','r');
formatSpec = '%f %f';
N_POM = fscanf(fileID,formatSpec,  [2 last_value])';
y_N_POM = N_POM(startValue:end,2)/average_factor;
x_N_POM = N_POM(startValue:end,1);

if(isNNZ)
fileID = fopen( root + 'N_POMconcVectorNNZ.txt','r');
formatSpec = '%f %f';
N_POM = fscanf(fileID,formatSpec,  [2 last_value])';
y_N_POMNNZ = N_POM(startValue:end,2);
x_N_POMNNZ = N_POM(startValue:end,1);
end


fileID = fopen( root + 'leakedNVector.txt','r');
formatSpec = '%f %f';
leakedN = fscanf(fileID,formatSpec,  [2 Inf])';
y_sumleakedC_N = leakedN(startValue:end,2)/average_factor;
x_leakedN = leakedN(startValue:end,1);

if(isNNZ)
fileID = fopen( root + 'leakedNVectorNNZ.txt','r');
formatSpec = '%f %f';
leakedNNNZ = fscanf(fileID,formatSpec,  [2 last_value])';
y_leakedNNNZ = leakedNNNZ(startValue:end,2);
x_leakedNNNZ = leakedNNNZ(startValue:end,1);
end

x = x_N_B;
plot(x,y_N_B./soilfactor,markers(file),'LineWidth', linewidth, 'Color',colors.microbial_biomass)
hold on
% plot(x,y_N_MN_old./soilfactor,'LineWidth', linewidth)
% hold on
% plot(x,y_N_MN_new./soilfactor,'LineWidth', linewidth)
% hold on
plot(x,y_N_POM./soilfactor,markers(file),'LineWidth', linewidth, 'Color',colors.pom)
hold on
plot(x,y_N_S./soilfactor,markers(file),'LineWidth', linewidth, 'Color',colors.dissolved_substrate)
hold on
plot(x,y_sumleakedC_N./soilfactor,markers(file),'LineWidth', linewidth, 'Color',colors.leakage_co2)
hold on
plot(x,(y_N_MN_new + y_N_MN_old)./soilfactor,markers(file),'LineWidth', linewidth, 'Color',colors.sum_microbial_necromass)
hold on
qw{file} = plot(x,x* 0,markers(file),'LineWidth', 0.01,'Color',[0 0 0]);
name = {'biomass' + names_appendix(file),  ...
    'POM'+names_appendix(file), 'dissolved substrate'+names_appendix(file),'leaked N'+names_appendix(file), ...
    'necromass total'+names_appendix(file)};
names = [names,name];

end %files
%names = {'biomass', 'Necromass old','Necromass new', 'POM', 'dissolved substrate', 'leaked N', 'necromass total'};
%legend(names, 'FontSize', 14, 'Location','southeast')
legend([qw{:}],  names_appendix, 'FontSize', fontsizes.legend, 'Location','southeast')
ax = gca
ax.XAxis.FontSize = fontsizes.xaxis
ax.YAxis.FontSize = fontsizes.yaxis
xlabel('days', 'FontSize', fontsizes.xlabel)

ylabel('amount nitrogen mg N g soil^{-1} (log scale)', 'FontSize', fontsizes.ylabel)

set(gca, 'YScale', 'log');  % Set y-axis to log scale
set(gca,'TickDir',tick_dir);
%set(gca, 'FontSize', 14)
if title_on
    title(titel, 'Interpreter', 'none', 'FontSize',fontsizes.title)
end
set(figure6, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure6, fullfile(outputfolder, "N_log_" + description + ".png"));
savefig(fullfile(outputfolder, "N_log_" + description))






%% flux data


figure10 = figure('visible', isVisible) 
names = {};
for file = 1:numel(files)
    appendix = files(file);
    %root = "/home.local/roetzer/C_N/txtdata_" + appendix +"/"
    root = "/home.other/fauamlx5g/roetzer/C_N/txtdata_" + appendix +"/"
    

fileID = fopen( root + 'CUE.txt','r');
formatSpec = '%f %f';
CUE = fscanf(fileID,formatSpec,  [2 last_value])';
y_CUE = CUE(startValue:end,2)/average_factor;
x_CUE = CUE(startValue:end,1);

fileID = fopen( root + 'f_B_C.txt','r');
formatSpec = '%f %f';
f_B_C = fscanf(fileID,formatSpec,  [2 last_value])';
y_f_B_C = f_B_C(startValue:end,2)/average_factor;
x_f_B_C = f_B_C(startValue:end,1);

fileID = fopen( root + 'f_B_N.txt','r');
formatSpec = '%f %f';
f_B_N = fscanf(fileID,formatSpec,  [2 last_value])';
y_f_B_N = f_B_N(startValue:end,2)/average_factor;
x_f_B_N = f_B_N(startValue:end,1);

fileID = fopen( root + 'R.txt','r');
formatSpec = '%f %f';
R = fscanf(fileID,formatSpec,  [2 last_value])';
y_R = R(startValue:end,2)/average_factor;
x_R = R(startValue:end,1);

fileID = fopen( root + 'R_O.txt','r');
formatSpec = '%f %f';
R_O = fscanf(fileID,formatSpec,  [2 last_value])';
y_R_O = R_O(startValue:end,2)/average_factor;
x_R_O = R_O(startValue:end,1);

fileID = fopen( root + 'f_BD_C.txt','r');
formatSpec = '%f %f';
f_BD_C = fscanf(fileID,formatSpec,  [2 last_value])';
y_f_BD_C = f_BD_C(startValue:end,2)/average_factor;
x_f_BD_C = f_BD_C(startValue:end,1);

fileID = fopen( root + 'f_BD_N.txt','r');
formatSpec = '%f %f';
f_BD_N = fscanf(fileID,formatSpec,  [2 last_value])';
y_f_BD_N = f_BD_N(startValue:end,2)/average_factor;
x_f_BD_N = f_BD_N(startValue:end,1);


fileID = fopen( root + 'f_C.txt','r');
formatSpec = '%f %f';
f_C = fscanf(fileID,formatSpec,  [2 last_value])';
y_f_C = f_C(startValue:end,2)/average_factor;
x_f_C = f_C(startValue:end,1);

fileID = fopen( root + 'f_N.txt','r');
formatSpec = '%f %f';
f_N = fscanf(fileID,formatSpec,  [2 last_value])';
y_f_N = f_N(startValue:end,2)/average_factor;
x_f_N = f_N(startValue:end,1);



fileID = fopen( root + 'f_MN_C.txt','r');
formatSpec = '%f %f';
f_MN_C = fscanf(fileID,formatSpec,  [2 last_value])';
y_f_MN_C = f_MN_C(startValue:end,2)/average_factor;
x_f_MN_C = f_MN_C(startValue:end,1);

fileID = fopen( root + 'f_MN_N.txt','r');
formatSpec = '%f %f';
f_MN_N = fscanf(fileID,formatSpec,  [2 last_value])';
y_f_MN_N = f_MN_N(startValue:end,2)/average_factor;
x_f_MN_N = f_MN_N(startValue:end,1);



fileID = fopen( root + 'f_POM_C.txt','r');
formatSpec = '%f %f';
f_POM_C = fscanf(fileID,formatSpec,  [2 last_value])';
y_f_POM_C = f_POM_C(startValue:end,2)/average_factor;
x_f_POM_C = f_POM_C(startValue:end,1);



x = x_f_B_C;

plot(x_CUE(195:300),y_R(195:300)./soilfactor,markers(file),'LineWidth', linewidth, 'Color',colors.r)
hold on 
plot(x_CUE(195:300),y_R_O(195:300)./soilfactor,markers(file),'LineWidth', linewidth, 'Color',colors.r_o)
%plot(x_CUE(195:300),y_R_O(195:300)/y_R(195:300),markers(file),'LineWidth', linewidth, 'Color',colorCode(2,:))
hold on 
% plot(x_CUE,y_f_C./soilfactor,'LineWidth', linewidth)
% hold on 
% plot(x_CUE,y_f_BD_C./soilfactor,'LineWidth', linewidth)
% hold on 
% plot(x_CUE,y_f_MN_C./soilfactor,'LineWidth', linewidth)
% hold on 
% plot(x_CUE,y_f_POM_C./soilfactor, markers(file),'LineWidth', linewidth, 'Color',colorCode(2,:))
% hold on
qw{file} = plot(x(195:300),(y_C_MN_new(195:300) + y_C_MN_old(195:300))./soilfactor * 0,markers(file),'LineWidth', 0.01,'Color',[0 0 0]);
hold on

name = {"R"+names_appendix(file),"R_O"+names_appendix(file),...
    % "f_C"+names_appendix(file),"f_BD_C"+names_appendix(file),...
    % "f_Nec_C"+names_appendix(file), "f_POM_C"+names_appendix(file)};
    };
names = [names,name];
f_POM_increae = max(y_f_POM_C)/mean(y_f_POM_C(100:199))
end %files
%legend(names, 'FontSize', 14, 'Location','southeast', 'Interpreter', 'none')
legend([qw{:}],  names_appendix, 'FontSize', fontsizes.legend, 'Location','southeast')
ax = gca
ax.XAxis.FontSize = fontsizes.xaxis
ax.YAxis.FontSize = fontsizes.yaxis

xlabel('days', 'FontSize', fontsizes.xlabel)
ylabel('amount carbon mg C g soil^{-1} per day ', 'FontSize', fontsizes.ylabel)
%set(gca, 'FontSize', 14)
set(gca,'TickDir',tick_dir);
if title_on
    title(titel, 'Interpreter', 'none', 'FontSize',fontsizes.title)
end
set(figure10, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure10, fullfile(outputfolder, "fluxesC_" + description + ".png"));
savefig(fullfile(outputfolder, "fluxesC_" + description))

figure11 = figure('visible', isVisible) 
names = {};
for file = 1:numel(files)
appendix = files(file);
%root = "/home.local/roetzer/C_N/txtdata_" + appendix +"/"
root = "/home.other/fauamlx5g/roetzer/C_N/txtdata_" + appendix +"/"

fileID = fopen( root + 'f_B_N.txt','r');
formatSpec = '%f %f';
f_B_N = fscanf(fileID,formatSpec,  [2 last_value])';
y_f_B_N = f_B_N(startValue:end,2)/average_factor;
x_f_B_N = f_B_N(startValue:end,1);


fileID = fopen( root + 'f_BD_N.txt','r');
formatSpec = '%f %f';
f_BD_N = fscanf(fileID,formatSpec,  [2 last_value])';
y_f_BD_N = f_BD_N(startValue:end,2)/average_factor;
x_f_BD_N = f_BD_N(startValue:end,1);


fileID = fopen( root + 'f_N.txt','r');
formatSpec = '%f %f';
f_N = fscanf(fileID,formatSpec,  [2 last_value])';
y_f_N = f_N(startValue:end,2)/average_factor;
x_f_N = f_N(startValue:end,1);



fileID = fopen( root + 'f_MN_N.txt','r');
formatSpec = '%f %f';
f_MN_N = fscanf(fileID,formatSpec,  [2 last_value])';
y_f_MN_N = f_MN_N(startValue:end,2)/average_factor;
x_f_MN_N = f_MN_N(startValue:end,1);



fileID = fopen( root + 'f_POM_N.txt','r');
formatSpec = '%f %f';
f_POM_N = fscanf(fileID,formatSpec,  [2 last_value])';
y_f_POM_N = f_POM_N(startValue:end,2)/average_factor;
x_f_POM_N = f_POM_N(startValue:end,1);

fileID = fopen( root + 'sumleakedN_S.txt','r');
formatSpec = '%f %f';
R_leaked = fscanf(fileID,formatSpec,  [2 last_value])';
y_R_leaked = R_leaked(startValue:end,2)/average_factor;
x_R_leaked = R_leaked(startValue:end,1);


X = x_f_POM_N;
plot(x_CUE,y_R_leaked./soilfactor,markers(file),'LineWidth', linewidth, 'Color',colors.r)
hold on 
plot(x_CUE,y_f_N./soilfactor,markers(file),'LineWidth', linewidth, 'Color',colors.dissolved_substrate)
hold on 
plot(x_CUE,y_f_BD_N./soilfactor,markers(file),'LineWidth', linewidth, 'Color',colors.microbial_biomass)
hold on 
plot(x_CUE,y_f_MN_N./soilfactor,markers(file),'LineWidth', linewidth, 'Color',colors.new_microbial_necromass)
hold on 
plot(x_CUE,y_f_POM_N./soilfactor,markers(file),'LineWidth', linewidth, 'Color',colors.pom)
hold on
qw{file} = plot(x,(y_C_MN_new + y_C_MN_old)./soilfactor * 0,markers(file),'LineWidth', 0.01,'Color',[0 0 0]);

name = {"R_leaked"+names_appendix(file), "f_N"+names_appendix(file), ...
    "f_BD_N"+names_appendix(file), "f_Nec_N"+names_appendix(file), "f_POM_N"+names_appendix(file)
    };
names = [names,name];

f_POM_increase = max(y_f_POM_N)/mean(y_f_POM_N(100:199))
end
ax = gca
ax.XAxis.FontSize = fontsizes.xaxis
ax.YAxis.FontSize = fontsizes.yaxis
xlabel('days', 'FontSize', fontsizes.xlabel)
ylabel('amount nitrogen mg N g soil^{-1} per day', 'FontSize', fontsizes.ylabel)
set(gca,'TickDir',tick_dir);

%set(gca,'ylim', [-2 1]);
%legend(names, 'FontSize', 14, 'Location','southeast', 'Interpreter', 'none')
legend([qw{:}],  names_appendix, 'FontSize', fontsizes.legend, 'Location','southeast')
%set(gca, 'FontSize', 14)
if title_on
    title(titel, 'Interpreter', 'none', 'FontSize',fontsizes.title)
end
set(figure11, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure11, fullfile(outputfolder, "fluxesN_" + description + ".png"));
savefig(fullfile(outputfolder, "fluxesN_" + description))









figure8 = figure('visible', isVisible) 
names = {};
for file = 1:numel(files)
    appendix = files(file);
    %root = "/home.local/roetzer/C_N/txtdata_" + appendix +"/"
    root = "/home.other/fauamlx5g/roetzer/C_N/txtdata_" + appendix +"/"

    
fileID = fopen( root + 'R.txt','r');
formatSpec = '%f %f';
R = fscanf(fileID,formatSpec,  [2 last_value])';
y_R = R(startValue:end,2)/average_factor;
x_R = R(startValue:end,1);

fileID = fopen( root + 'f_C.txt','r');
formatSpec = '%f %f';
f_C = fscanf(fileID,formatSpec,  [2 last_value])';
y_f_C = f_C(startValue:end,2)/average_factor;
x_f_C = f_C(startValue:end,1);

x_R = x_f_C;

CUE_A = 1 - y_R./abs(y_f_C);
plot(x_CUE,CUE_A,markers(file),'LineWidth', linewidth, 'Color',colors.cue)
hold on 
mean(CUE_A(100:200))
name = {"CUE"+names_appendix(file)};
names = [names,name];


end
%CUE_B = 1 - (y_R+ y_f_BD_C)./abs(y_f_C);
%plot(x_CUE,CUE_B,'LineWidth', linewidth)
ax = gca;
ax.XAxis.FontSize = fontsizes.xaxis;
ax.YAxis.FontSize = fontsizes.yaxis;
xlabel('days', 'FontSize', fontsizes.xlabel)
ylabel('Carbon Use Efficiency Microbes', 'FontSize', fontsizes.ylabel)
set(gca,'TickDir',tick_dir);

legend(names, 'FontSize', fontsizes.legend, 'Location','southeast', 'Interpreter', 'none')
%set(gca,'ylim', [-2 1]);
%set(gca,'xlim', [0 numel(x_CUE)]);
%legend(["CUE_A", "CUE_B"], 'Interpreter', 'none')
%set(gca, 'FontSize', 14)
if title_on
    title(titel, 'Interpreter', 'none','FontSize',fontsizes.title)
end
set(figure8, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure8, fullfile(outputfolder, "CUE_" +description + ".png"));
savefig(fullfile(outputfolder, "CUE_" +description))
end

if true
%% figure 10
startValue = 1;
figure10 = figure('visible', isVisible) 
add_bar = [];
timesteps = [200,201, 500];
y_compare = [];
y_compare_t1 = [];
y_compare_t2 = [];
y_compare_t3 = [];
for file = 1:numel(files)
appendix = files(file);
%root = "/home.local/roetzer/C_N/txtdata_" + appendix +"/"
root = "/home.other/fauamlx5g/roetzer/C_N/txtdata_" + appendix +"/"
    
%% carbon data
fileID = fopen(root + 'C_BVector.txt','r');
formatSpec = '%f %f';
C_B = fscanf(fileID,formatSpec,  [2 last_value])';
y_C_B = C_B(startValue:end,2)/average_factor;
x_C_B = C_B(startValue:end,1);

fileID = fopen(root + 'C_SVector.txt','r');
formatSpec = '%f %f';
C_S = fscanf(fileID,formatSpec,  [2 last_value])';
y_C_S = C_S(startValue:end,2)/average_factor;
x_C_S = C_S(startValue:end,1);


fileID = fopen( root + 'C_MNVector.txt','r');
formatSpec = '%f %f';
C_MN = fscanf(fileID,formatSpec,  [2 last_value])';
y_C_MN = C_MN(startValue:end,2)/average_factor;
x_C_MN = C_MN(startValue:end,1);

fileID = fopen( root + 'C_MNVector_old.txt','r');
formatSpec = '%f %f';
C_MN_old = fscanf(fileID,formatSpec,  [2 last_value])';
y_C_MN_old = C_MN_old(startValue:end,2)/average_factor;
x_C_MN_old = C_MN_old(startValue:end,1);

fileID = fopen( root + 'C_MNVector_new.txt','r');
formatSpec = '%f %f';
C_MN_new = fscanf(fileID,formatSpec,  [2 last_value])';
y_C_MN_new = C_MN_new(startValue:end,2)/average_factor;
x_C_MN_new = C_MN_new(startValue:end,1);

fileID = fopen( root + 'C_POMconcVector.txt','r');
formatSpec = '%f %f';
C_POM = fscanf(fileID,formatSpec,  [2 last_value])';
y_C_POM = C_POM(startValue:end,2)/average_factor;
x_C_POM = C_POM(startValue:end,1);

fileID = fopen( root + 'CO2Vector.txt','r');
formatSpec = '%f %f';
CO2 = fscanf(fileID,formatSpec,  [2 last_value])';
y_CO2 = CO2(startValue:end,2)/average_factor;
x_CO2 = CO2(startValue:end,1);

fileID = fopen( root + 'CO2Vector_over.txt','r');
formatSpec = '%f %f';
CO2_over = fscanf(fileID,formatSpec,  [2 last_value])';
y_CO2_over = CO2_over(startValue:end,2)/average_factor;
x_CO2_over = CO2_over(startValue:end,1);

POM_added = floor([1:numel(y_C_B)]./10) * 0.083;

y_C = [];
y_C_sum = [];
CO2_days = [];
for i = 1:numel(y_C_MN)
   y_C_add = [y_C_B(i), y_C_MN_old(i),y_C_MN_new(i),y_C_POM(i),  y_C_S(i),y_CO2(i)];%
   sum_C = y_C_B(i)+ y_C_MN_old(i)+y_C_MN_new(i)+y_C_POM(i)+  y_C_S(i);
   y_C_add = y_C_add./soilfactor;% / soilParticleDensity;
   y_C = [y_C; y_C_add];
   y_C_sum = [y_C_sum; sum_C];

   if((i-1) == 0)
       CO2_day_ = 0;
       CO2_day = 0;
   else
      CO2_day_ = y_CO2(i) - y_CO2(i-1);
      CO2_day = CO2_day_./( soil_particleNNZ * soilParticleDensity/1000);
   end
   CO2_days = [CO2_days; CO2_day/24];
end
TOC_percent = y_C_sum./(y_C_sum + soil_particleNNZ * soilParticleDensity);
C_mic = y_C_B(i)./y_C_sum;
TOC = y_C_sum./(soil_particleNNZ * soilParticleDensity/1000); 
CO2_hour_avag = mean(CO2_days);
Biomass = y_C_B./soilfactor;

if(file == 1)
    y_compare = [y_compare;y_C(1,:) ]
    y_compare_t1 = [y_compare_t1; y_C(timesteps(1),:)];
    y_compare_t2 = [y_compare_t2; y_C(timesteps(2),:)];
end
    
    %y_compare_t2 = [y_compare_t2; y_C(timesteps(2),:)];
    y_compare_t3 = [y_compare_t3; y_C(timesteps(3),:)];
    %y_compare = [y_compare;y_C(199,:) ]


end % files
add_null = y_C(500,:) * 0;
y_compare = [y_compare;add_null; y_compare_t1; add_null; y_compare_t2;add_null; y_compare_t3]
add_bar = [add_bar,[2,4,6]];


% names_appendix = [names_appendix(1) + ' initial', "add POM1"...
%     ,names_appendix(1) + ' day199', "add exudation"...
%     ,names_appendix(1) + ' day200', "add POM2"...
%     ,names_appendix + ' day500']

names_xlabel = ['day 0', sprintf('addition\\newlineday 1-199')...
    ,'day 199', sprintf('addition\\newlineday 200')...
    ,'day 200', sprintf('addition\\newlineday 201 -500')...
    ,names_appendix + sprintf('\\newlineday 500')]

% names_xlabel = {'day 0', sprintf('addition\\newlineday 1-199')...
%     ,'day 199', sprintf('addition\\newlineday 200')...
%     ,'day 200', sprintf('addition\\newlineday 201 -500')...
%     ,sprintf('%s\\newlineday 500',names_appendix )}

b1 = bar(y_compare, 'stacked');
names = {'biomass',  'Necromass old', 'Necromass new', 'POM', 'dissolved substrate','CO2'};
colornames = {'microbial_biomass',  'old_microbial_necromass', 'new_microbial_necromass', 'pom', 'dissolved_substrate','leakage_co2'};
opt_ax = [0.1300 0.1476 0.7750 0.7774];

ax = gca;
ax.Units ='normalized';
ax.Position = opt_ax;
ax.ActivePositionProperty = 'position';
ax.XAxis.FontSize = fontsizes.xaxis;
ax.XTickLabel = names_xlabel
ax.YAxis.FontSize = fontsizes.yaxis;
set(figure10, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

for i = 1:numel(b1)
    b1(i).FaceColor = 'flat';
    b1(i).CData = getfield(colors, colornames{i});
end

xlabel('days', 'FontSize', fontsizes.xlabel)
%ylabel('amount carbon g cm^{-3}', 'FontSize', 14)
ylabel('amount carbon mg C g soil^{-1}', 'FontSize', fontsizes.ylabel)
%set(gca,'ylim',ylim_C_stacked);
if title_on
    title(titel , 'Interpreter', 'none','FontSize',fontsizes.title)
end
%set(gca, 'FontSize', 14)
percentage_plot = 1;
% Add percentage annotations
for k = 1:size(y_compare, 1)
    for j = 1:size(y_compare, 2)
        if mod(k, percentage_plot) == 0 && ~ismember(k,add_bar) 
            percentage = y_compare(k, j) / sum(y_compare(k, :)) * 100;
            if y_compare(k, j) ==0
                continue
            end
            tmp_color = 'black';
            if y_compare(k, j) <0.03 && percentage ~=0 
                

                
                End   = [k, sum(y_compare(k, 1:j)) - y_compare(k, j)/2]  % x , y
                Start = [End(1)-0.15, End(2)+0.5]; % x, y
                Offset = 0.35

                xapf = @(x,pos,xl) pos(3)*(x-(min(xl)))/(diff(xl))+pos(1);                % 'x' Annotation Position Function
                yapf = @(y,pos,yl) pos(4)*(y-min(yl))/diff(yl)+pos(2);                % 'y' Annotation Position Function
                xl = [0 length(names_xlabel)+1+Offset]%double(xlim) + [-1 1]

                
                yl = ylim
                %pause
                pos = opt_ax
                ax
                %annotation('textarrow', [0.1 0.1300], [0.13 0.1476])
                %a=0.1476+0.7774
                %annotation('textarrow', [0.1 0.1300], [0.13 a])
                annotation('textarrow', xapf([Start(1) End(1)+Offset/2],pos,xl), yapf([Start(2) End(2)],pos,yl))
                xVal=1
                yVal=5
                
                
                set(figure10, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
                text(Start(1)-Offset, Start(2), sprintf('%.2f\n(%i%%)',y_compare(k, j), ceil(percentage)), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', fontsizes.bartext_arrow,'Color', tmp_color)%12
                continue
            end
            if strcmp(colornames{j}, 'leakage_co2') || strcmp(colornames{j}, 'dissolved_substrate')
                tmp_color = colors.leakage_co2_text;
            end
            text(k, sum(y_compare(k, 1:j)) - y_compare(k, j)/2, sprintf('%.2f (%i%%)',y_compare(k, j), ceil(percentage)), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', fontsizes.bartext,'Color', tmp_color)%12
        end
    end
end

hold on 
size_overlay = size(y_compare);
y_add = zeros([size_overlay(1), 3]);
y_add(add_bar(1),1) = sum(y_compare(add_bar(1)-1,:));
y_add(add_bar(1),2) = POM_added(timesteps(1)-1);

y_add(add_bar(2),1) = sum(y_compare(add_bar(2)-1,:));%y_C(timesteps(1)
y_add(add_bar(2),2) = POM_added(timesteps(2)) - POM_added(timesteps(1)-1);
y_add(add_bar(2),3) = (sum(y_compare(add_bar(3)-1,:)) - sum(y_compare(add_bar(2)-1,:))) -(POM_added(timesteps(2)) - POM_added(timesteps(1)-1));
%sum(y_C(timesteps(2),:)) - sum(y_C(timesteps(1),:))
y_add(add_bar(3),1) = sum(y_compare(add_bar(3)-1,:));%y_C(timesteps(2)
y_add(add_bar(3),2) = POM_added(499) - POM_added(200);
y_add(add_bar(3),3) = (sum(y_compare(end,:)) - sum(y_compare(end-1,:)));




b2 = bar(y_add, 'stacked')%, 'FaceColor','flat'

% make hatched lines
epsX = 0.005;
epsY = 0.05;
barWidth = 0.4-epsX; % Balkenbreite fÃ¼r korrektes Platzieren der Linien
for j = 2:3
    for i = 1:length(b2(j).XData)  % Iteriere Ã¼ber die Balken
        xPos = i + eps;  % x-Position des Balkens
        yHeight = y_add(i,j);   % HÃ¶he des ersten Balkensegments
        yStartHeight=0;
        for k = j-1:-1:1
            yStartHeight = yStartHeight + y_add(i,k) + epsY;
        end
        % Schraffur-Linien innerhalb des Balkens zeichnen
        for yStripe = yStartHeight:0.25:(yHeight+yStartHeight-2*epsY)
    
            x1 = xPos - barWidth;
            x2 = xPos + barWidth;
            y1 = yStripe;
            y2 = yStripe;
            orgColor = 'black';
            alpha = 0.5;
            switch j
                case 2
                    orgColor = getfield(colors, colornames{4});
                case 3
                    orgColor = getfield(colors, colornames{5});
            end
            hatchColor = (1-alpha) * orgColor + alpha * [1 1 1];
            plot([x1, x2], ...
                [y1, y2], 'k', 'LineWidth', 2, 'Color',hatchColor);
        end
    end
end


% b2(1).AlphaData(2,:) = [1 1 1];
% b2(2).CData(2,:) = colorCode(4,:);
% b2(3).CData(2,:) = colorCode(5,:);
b2(1).FaceColor = 'none';%[1 1 1];
b2(1).EdgeColor = [1 1 1];
%b2(2).FaceColor = colorCode(4,:);
b2(2).FaceColor = getfield(colors, colornames{4});
b2(2).EdgeColor = 'black';
%b2(3).FaceColor = colorCode(5,:);
b2(3).FaceColor = getfield(colors, colornames{5});
b2(3).EdgeColor = 'black';



for k = 1:size(y_add, 1)
    for j = 1:size(y_add, 2)
       if(y_add(k, j) > 0 &&j  > 1)
        text(k, sum(y_add(k, 1:j)) - y_add(k, j)/2, sprintf('%.2f \n ',y_add(k, j)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', fontsizes.bartext)%12
        end
    end
end

%set(figure10, 'Units', 'normalized', 'OuterPosition', [0 0 1 1], 'Position', opt_ax);
%[1 0.1841 0.1841]
legend(names, 'FontSize', fontsizes.legend, 'Location','northwest')
set(gca,'TickDir',tick_dir);
saveas(figure10, fullfile(outputfolder, "C_initial_final_" + description + ".png"));
savefig(fullfile(outputfolder, "C_initial_final_" + description))
ax = gca().get
%pause;
end


















system("mogrify -trim " + outputfolder + "*.png")




close all











end
% fileID = fopen( root + 'CUE.txt','r');
% formatSpec = '%f %f';
% CUE = fscanf(fileID,formatSpec,  [2 last_value])';
% y_CUE = CUE(startValue:end,2)/average_factor;
% x_CUE = CUE(startValue:end,1);




% figure13= figure('visible', isVisible) 
% C_N = y_C_S./y_N_S;
% plot(x_CUE,C_N,'LineWidth', linewidth)
% xlabel('days', 'FontSize', 14)
% ylabel('C_N ', 'FontSize', 14)
% legend(["C_N"], 'Interpreter', 'none')
% set(gca, 'FontSize', 14)
% title(titel, 'Interpreter', 'none')
% set(figure13, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
% saveas(figure13, fullfile(outputfolder, "C_N_" + appendix + ".png"));