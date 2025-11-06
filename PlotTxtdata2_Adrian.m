function PlotTxtdata2_Adrian(file, description)
close all
if nargin < 2
    description =file;
end
if nargin < 1
    file = "noMove_longterm";
end
imageFolder = "bilder/CN/paperReady/";
%imageFolder = "final_pics2/";
appendix = file%"noMove_longterm_noreactiveEdgesforDecay" %"noMove_longterm" %mucilageC_after100stepsCN01
outputfolder = imageFolder + appendix +"/";
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
%root = "/home.local/roetzer/C_N/txtdata_" + appendix +"/"
root = "/home.other/fauamlx5g/roetzer/C_N/txtdata_" + appendix +"/"
titel = description;
percentage_plot = 50;%;%25;
isVisible = 'on';
average_factor = 1;
if file.contains("2years")
    last_value = 2000;
else
last_value = Inf;
end
% DOC ist zu hoch

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
    'xlabel', 26, ...
    'ylabel', 26, ...
    'xaxis', 26, ...
    'yaxis', 26, ... 
    'legend', 26, ...
    'title', 30, ...
    'bartext', 18, ...
    'bartext_arrow', 18 ...
);
% (days mod mod_days) for text in stacked bar graph
mod_days = 50;

title_on = false;
tick_dir = 'out'; %or 'in'


ylim_C_stacked = [0 0.2];
ylim_N_stacked = [0 4*10^-3];

ylim_C = [10^-8 10^-0];
ylim_N = [10^-8 10^0 ];
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
y_C = [];
y_C_sum = [];
CO2_days = [];
for i = 1:numel(y_C_MN)
   y_C_add = [y_C_B(i), y_C_MN_old(i),y_C_MN_new(i),y_C_POM(i),  y_C_S(i)];%,y_CO2(i)
   %y_C_add = [y_C_B(i)/y_C_BNNZ(i), y_C_S(i)/y_C_SNNZ(i),y_C_MN_old(i)/y_C_MN_oldNNZ(i),y_C_MN_new(i)/y_C_MN_newNNZ(i), y_C_POM(i)/y_C_POMNNZ(i), y_CO2(i)/y_CO2NNZ(i)];
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


x = x_N_MN;
y_N = [];
y_N_sum = [];
for i = 1:numel(y_N_MN)
   %y_add = [y_N_B(i), y_N_S(i), y_N_MN(i), y_N_POM(i), y_sumleakedC_N(i)];
   y_N_add = [y_N_B(i), y_N_MN_old(i),y_N_MN_new(i), y_N_POM(i),y_N_S(i)];%,y_sumleakedC_N(i)
   %y_N_add = [y_N_B(i)/y_N_BNNZ(i), y_N_MN_old(i)/y_N_MN_oldNNZ(i),y_N_MN_new(i)/y_N_MN_newNNZ(i), y_N_POM(i)/y_N_POMNNZ(i),y_N_S(i)/y_N_SNNZ(i),  y_leakedN(i)/y_leakedNNNZ(i)];
   sum_N = y_N_B(i) + y_N_S(i) + y_N_MN(i) + y_N_POM(i);
   y_N_add = y_N_add./soilfactor;%/ soilParticleDensity;
   y_N = [y_N; y_N_add];
   y_N_sum = [y_N_sum; sum_N];
end
TON_percent = y_N_sum./(y_N_sum + soil_particleNNZ * soilParticleDensity);
TON = y_N_sum./(soil_particleNNZ * soilParticleDensity/1000); 
%% Plots 
percentage_plot = round(numel(x_N_POM)/20);%;%25;
linewidth = 4;
if true

%% figure 1
figure1 = figure('visible', isVisible)
names = {'biomass', 'dissolved substrate', 'Necromass', 'POM', 'CO2'};
names = {'biomass',  'Necromass old', 'Necromass new', 'POM', 'dissolved substrate'};%,'CO2'
colornames = {'microbial_biomass',  'old_microbial_necromass', 'new_microbial_necromass', 'pom', 'dissolved_substrate'};

b1 = bar(x(1:end), y_C, 'stacked')
for i = 1:numel(b1)
    b1(i).FaceColor = 'flat';
    b1(i).CData = getfield(colors, colornames{i})
end
legend(names, 'FontSize', fontsizes.legend, 'Location','northwest')
ax = gca;
ax.XAxis.FontSize = fontsizes.xaxis;
ax.YAxis.FontSize = fontsizes.yaxis;
xlabel('days', 'FontSize', fontsizes.xlabel)
%ylabel('amount carbon g cm^{-3}', 'FontSize', 14)
ylabel('amount carbon mg C g soil^{-1}', 'FontSize', fontsizes.ylabel)
%set(gca,'ylim',ylim_C_stacked);
if title_on
    title(titel , 'Interpreter', 'none','FontSize',fontsizes.title)
end
%set(gca, 'FontSize', 14)

% Add percentage annotations
temp_counter = 0;
for k = 1:size(y_C, 1)
    temp_counter = temp_counter +1;
    for j = 1:size(y_C, 2)
        if mod(temp_counter, mod_days) == 0 %|| mod(k, percentage_plot) == 0 
            percentage = y_C(k, j) / sum(y_C(k, :)) * 100;
            if y_C(k, j) ==0
                continue
            end
            tmp_color = 'black';
            if strcmp(colornames{j}, 'leakage_co2')
                tmp_color = colors.leakage_co2_text;
            end
            vertAlign = 'middle';
            if percentage <= 3
                vertAlign = 'bottom';
            end
            text(x(k), sum(y_C(k, 1:j)) - y_C(k, j)/2, sprintf('%.2f\n(%i%%)',y_C(k, j), ceil(percentage)), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', vertAlign, 'FontSize', fontsizes.bartext,'Color', tmp_color)%12

        end
    end
end
set(gca,'TickDir',tick_dir);
set(figure1, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure1, fullfile(outputfolder, "C_dist_" + appendix + ".png"));
savefig(fullfile(outputfolder, "C_dist_" + appendix))


%% figure 2
figure2 = figure('visible', isVisible)
names = {'biomass', 'dissolved substrate', 'Necromass', 'POM', 'leaked N'};
names = {'biomass', 'Necromass old','Necromass new', 'POM', 'dissolved substrate'};%, 'leaked N'
colornames = {'microbial_biomass',  'old_microbial_necromass', 'new_microbial_necromass', 'pom', 'dissolved_substrate'};

b2 = bar(x(1:end), y_N, 'stacked')
for i = 1:numel(b2)
    b2(i).FaceColor = 'flat';
    b2(i).CData = getfield(colors, colornames{i})
end
legend(names, 'FontSize', fontsizes.legend, 'Location','northwest')
ax = gca;
ax.XAxis.FontSize = fontsizes.xaxis;
ax.YAxis.FontSize = fontsizes.yaxis;
xlabel('days', 'FontSize', fontsizes.xlabel)
%ylabel('amount nitrogen g cm^{-3}', 'FontSize', 14)
ylabel('amount nitrogen mg N g soil^{-1}', 'FontSize', fontsizes.ylabel)
if title_on
    title(titel, 'Interpreter', 'none','FontSize', fontsizes.title)
end
%set(gca,'ylim',ylim_N_stacked);
%set(gca, 'FontSize', 14)

% Add percentage annotations
temp_counter = 0;
for k = 1:size(y_N, 1)
    temp_counter = temp_counter + 1;
    for j = 1:size(y_N, 2)
        if mod(temp_counter, mod_days) == 0 % if mod(k, percentage_plot) == 0 || k == 1
            percentage = y_N(k, j) / sum(y_N(k, :)) * 100;
            vertAlign = 'middle';
            if percentage <= 3
                vertAlign = 'bottom';
            end
            text(x(k), sum(y_N(k, 1:j)) - y_N(k, j)/2, sprintf('%.2f\n(%i%%)',y_N(k, j), ceil(percentage)), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', vertAlign, 'FontSize', fontsizes.bartext)
        end
    end
end
set(gca,'TickDir',tick_dir);
set(figure2, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure2, fullfile(outputfolder, "N_dist_" + appendix + ".png"));
savefig(fullfile(outputfolder, "N_dist_" + appendix))
%% figure 5
% smallest_value = 10^-2;
% threshold = 10^-2;
% y_N_B(find(y_N_B <=threshold,1)) = smallest_value;
% y_N_MN_old(find(y_N_MN_old <=threshold,1))= smallest_value;
% y_N_MN_new(find(y_N_MN_new <=threshold,1))= smallest_value;
% y_N_POM(find(y_N_POM <=threshold,1))= smallest_value;
% y_N_S( find(y_N_S<=threshold,1)) = smallest_value;
% y_sumleakedC_N(find(y_sumleakedC_N <=threshold,1))= smallest_value;
% 
% 
% y_N_B(y_N_B <threshold) = NaN;
% y_N_MN_old(y_N_MN_old <threshold)= NaN;
% y_N_MN_new(y_N_MN_new <threshold)= NaN;
% y_N_POM(y_N_POM <threshold)= NaN;
% y_N_S( y_N_S<threshold) = NaN;
% y_sumleakedC_N(y_sumleakedC_N <threshold)= NaN;

figure5 = figure('visible', isVisible) 
names = {'biomass', 'Necromass old','Necromass new', 'POM', 'dissolved substrate', 'leaked N', 'necromass total'};
linewidth = 4;
plot(x,y_N_B./soilfactor,'LineWidth', linewidth, 'Color', colors.microbial_biomass)
hold on
plot(x,y_N_MN_old./soilfactor,'LineWidth', linewidth, 'Color', colors.old_microbial_necromass)
hold on
plot(x,y_N_MN_new./soilfactor,'LineWidth', linewidth, 'Color', colors.new_microbial_necromass)
hold on
plot(x,y_N_POM./soilfactor,'LineWidth', linewidth, 'Color', colors.pom)
hold on
plot(x,y_N_S./soilfactor,'LineWidth', linewidth, 'Color', colors.dissolved_substrate)
hold on
plot(x,y_sumleakedC_N./soilfactor,'LineWidth', linewidth, 'Color', colors.leakage_co2)
hold on
plot(x,(y_N_MN_new + y_N_MN_old)./soilfactor,'LineWidth', linewidth, 'Color', colors.sum_microbial_necromass)

legend(names, 'FontSize', fontsizes.legend, 'Location','southeast')
ax = gca;
ax.XAxis.FontSize = fontsizes.xaxis;
ax.YAxis.FontSize = fontsizes.yaxis;
xlabel('days', 'FontSize', fontsizes.xlabel)
%ylabel('amount nitrogen g cm^{-3} (log scale)', 'FontSize', 14)
ylabel('amount nitrogen mg N g soil^{-1} ', 'FontSize', fontsizes.ylabel)%(log scale)
% if(~(file.contains("normal") || file.contains("noPOMDecay")))
% set(gca,'ylim',ylim_N);
% end
%set(gca, 'YScale', 'log');  % Set y-axis to log scale


%ylim = mylimits;
%set(gca, 'FontSize', 14)
if title_on
    title(titel, 'Interpreter', 'none','FontSize',fontsizes.title)
end
set(gca,'TickDir',tick_dir);
set(figure5, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure5, fullfile(outputfolder, "N_log_" + appendix + ".png"));
savefig(fullfile(outputfolder, "N_log_" + appendix))

%% figure 6
figure6 = figure('visible', isVisible) 
linewidth = 4;

% smallest_value = 10^-2;
% threshold = 10^-2;
% y_C_B(find(y_C_B <=threshold,1)) = smallest_value;
% y_C_MN_old(find(y_C_MN_old <=threshold,1))= smallest_value;
% y_C_MN_new(find(y_C_MN_new <=threshold,1))= smallest_value;
% y_C_POM(find(y_C_POM <=threshold,1))= smallest_value;
% y_C_S( find(y_C_S<=threshold,1)) = smallest_value;
% y_CO2(find(y_CO2 <=threshold,1))= smallest_value;
% 
% 
% y_C_B(y_C_B <threshold) = NaN;
% y_C_MN_old(y_C_MN_old <threshold)= NaN;
% y_C_MN_new(y_C_MN_new <threshold)= NaN;
% y_C_POM(y_C_POM <threshold)= NaN;
% y_C_S( y_C_S<threshold) = NaN;
% y_CO2(y_CO2 <threshold)= NaN;

POM_added = floor([1:numel(y_C_B)]./10) * 0.083;

plot(x,y_C_B./soilfactor,'LineWidth', linewidth, 'Color', colors.microbial_biomass)
hold on
plot(x,y_C_MN_old./soilfactor,'LineWidth', linewidth, 'Color', colors.old_microbial_necromass)
hold on
plot(x,y_C_MN_new./soilfactor,'LineWidth', linewidth, 'Color', colors.new_microbial_necromass)
hold on
plot(x,y_C_POM./soilfactor,'LineWidth', linewidth, 'Color', colors.pom)
hold on
plot(x,y_C_S./soilfactor,'LineWidth', linewidth, 'Color', colors.dissolved_substrate)
hold on
plot(x,y_CO2./soilfactor,'LineWidth', linewidth, 'Color', colors.leakage_co2)
hold on
plot(x,(y_C_MN_new + y_C_MN_old)./soilfactor,'LineWidth', linewidth, 'Color', colors.sum_microbial_necromass)
hold on
%plot(x,POM_added,'LineWidth', linewidth, 'Color', colors.pom_added)
%hold on
%plot(x,y_CO2_over,'LineWidth', linewidth)

%names = {'biomass',  'Necromass old', 'Necromass new', 'POM', 'dissolved substrate','CO2', 'CO2_over'};
names = {'biomass',  'Necromass old', 'Necromass new', 'POM', 'dissolved substrate','CO2', 'necromass total'};%, 'POM added'
legend(names, 'FontSize', fontsizes.legend, 'Location','southeast')
ax = gca;
ax.XAxis.FontSize = fontsizes.xaxis;
ax.YAxis.FontSize = fontsizes.yaxis;
xlabel('days', 'FontSize', fontsizes.xlabel)
%ylim(mylimits);
%ylabel('amount carbon g cm^{-3} (log scale)', 'FontSize', 14)
ylabel('amount carbon mg C g soil^{-1} ', 'FontSize', fontsizes.ylabel)%(log scale)
% if(~(file.contains("normal") || file.contains("noPOMDecay")))
% set(gca,'ylim',ylim_C);
% end
%set(gca, 'YScale', 'log');  % Set y-axis to log scale
%ylim([-0.1,1.1])
%set(gca, 'FontSize', 14)
if title_on
    title(titel, 'Interpreter', 'none','FontSize',fontsizes.title)
end
set(gca,'TickDir',tick_dir);
set(figure6, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure6, fullfile(outputfolder, "C_log_" + appendix + ".png"));
%print(fullfile(outputfolder, 'figure_large_dpi.png'), '-dpng', '-r300');
savefig(fullfile(outputfolder, "C_log_" + appendix))
%% figure 7
% figure7 = figure('visible', isVisible) 
% linewidth = 4;
% 
% plot(x,gradient(y_C_B./soilfactor),'LineWidth', linewidth)
% hold on
% plot(x,gradient(y_C_MN_old./soilfactor),'LineWidth', linewidth)
% hold on
% plot(x,gradient(y_C_MN_new./soilfactor),'LineWidth', linewidth)
% hold on
% plot(x,gradient(y_C_POM./soilfactor),'LineWidth', linewidth)
% hold on
% plot(x,gradient(y_C_S./soilfactor),'LineWidth', linewidth)
% hold on
% plot(x,gradient(y_CO2./soilfactor),'LineWidth', linewidth)
% hold on
% plot(x,gradient((y_C_MN_new + y_C_MN_old)./soilfactor),'LineWidth', linewidth)
% 
% names = {'biomass',  'Necromass old', 'Necromass new', 'POM', 'dissolved substrate','CO2', 'necromass total'};
% legend(names, 'FontSize', 14, 'Location','southeast')
% xlabel('days', 'FontSize', 14)
% ylabel('amount carbon mg C g soil^{-1} (log scale)', 'FontSize', 14)
% %set(gca, 'YScale', 'log');
% set(gca, 'FontSize', 14)
% title(titel, 'Interpreter', 'none')
% set(figure7, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
% saveas(figure7, fullfile(outputfolder, "C_log_gradient" + appendix + ".png"));


end

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



figure8 = figure('visible', isVisible) 
%plot(x_CUE,y_CUE,'LineWidth', linewidth)
%hold on 
CUE_A = 1 - y_R./abs(y_f_C);
plot(x_CUE,CUE_A,'LineWidth', linewidth, 'Color',colors.cue)
%hold on 
mean(CUE_A(100:200))
CUE_B = 1 - (y_R+ y_f_BD_C)./abs(y_f_C);
ax = gca;
ax.XAxis.FontSize = fontsizes.xaxis;
ax.YAxis.FontSize = fontsizes.yaxis;
%plot(x_CUE,CUE_B,'LineWidth', linewidth)
xlabel('days', 'FontSize', fontsizes.xlabel)
ylabel('Carbon Use Efficiency Microbes', 'FontSize', fontsizes.ylabel)
%set(gca,'ylim', [-2 1]);
%set(gca,'xlim', [0 numel(x_CUE)]);
%legend(["CUE_A", "CUE_B"], 'Interpreter', 'none')
%set(gca, 'FontSize', 14)
if title_on
    title(titel, 'Interpreter', 'none','FontSize',fontsizes.title)
end
set(gca,'TickDir',tick_dir);
set(figure8, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure8, fullfile(outputfolder, "CUE_" + appendix + ".png"));
savefig(fullfile(outputfolder, "CUE_" + appendix))



if false

figure9 = figure('visible', isVisible) 
%plot(x_CUE,y_CUE,'LineWidth', linewidth)
%hold on 
ax = gca;
ax.XAxis.FontSize = fontsizes.xaxis;
ax.YAxis.FontSize = fontsizes.yaxis;
xlabel('days', 'FontSize', fontsizes.xlabel)
plot(x_CUE,y_R./soilfactor,'LineWidth', linewidth, 'Color',colors.r)
%hold on 
%plot(x_CUE,y_R_O./soilfactor,'LineWidth', linewidth)
hold on 
plot(x_CUE,y_f_C./soilfactor,'LineWidth', linewidth, 'Color',colors.dissolved_substrate)
hold on 
plot(x_CUE,y_f_N./soilfactor,'LineWidth', linewidth, 'Color',colors.dissolved_substrate)
hold on 
plot(x_CUE,y_f_BD_C./soilfactor,'LineWidth', linewidth, 'Color',colors.microbial_biomass)
hold on 
plot(x_CUE,y_f_BD_N./soilfactor,'LineWidth', linewidth, 'Color',colors.microbial_biomass)
hold on 
plot(x_CUE,y_f_MN_C./soilfactor,'LineWidth', linewidth, 'Color',colors.new_microbial_necromass)
hold on 
plot(x_CUE,y_f_MN_N./soilfactor,'LineWidth', linewidth, 'Color',colors.new_microbial_necromass)
hold on 
plot(x_CUE,y_f_POM_C./soilfactor,'LineWidth', linewidth, 'Color',colors.pom)
hold on 
plot(x_CUE,y_f_POM_N./soilfactor,'LineWidth', linewidth, 'Color',colors.pom)

ylabel('g cm 3 ', 'FontSize', fontsizes.ylabel)



%set(gca,'ylim', [-2 1]);
legend(["R", "f_U", "f_N", "f_BD_C", "f_BD_N", "f_Nec_C", "f_MN_N", "f_POM_C", "f_POM_N"], 'Interpreter', 'none', 'FontSize', fontsizes.legend)
%set(gca, 'FontSize', 14)
if title_on
    title(titel, 'Interpreter', 'none','FontSize', fontsizes.title)
end
set(gca,'TickDir',tick_dir);
set(figure9, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure9, fullfile(outputfolder, "fluxes_" + appendix + ".png"));
savefig(fullfile(outputfolder, "fluxes_" + appendix))
end




figure10 = figure('visible', isVisible) 
%plot(x_CUE,y_CUE,'LineWidth', linewidth)
%hold on 

plot(x_CUE,y_R./soilfactor,'LineWidth', linewidth, 'Color',colors.r)
hold on 
y_R_O_abs =  y_R_O;
y_R_O_abs(y_R_O_abs > 0) = 0;
y_R_O_abs = abs(y_R_O_abs);
fract = y_R_O./y_R;
max(fract([100:200]))
plot(x_CUE,y_R_O./soilfactor,'LineWidth', linewidth, 'Color',colors.r_o)
hold on 
plot(x_CUE,y_f_C./soilfactor,'LineWidth', linewidth, 'Color',colors.dissolved_substrate)
hold on 
plot(x_CUE,y_f_BD_C./soilfactor,'LineWidth', linewidth, 'Color',colors.microbial_biomass)
hold on 
plot(x_CUE,y_f_MN_C./soilfactor,'LineWidth', linewidth, 'Color',colors.new_microbial_necromass)
hold on 
plot(x_CUE,y_f_POM_C./soilfactor,'LineWidth', linewidth, 'Color',colors.pom)
%hold on 
%plot(x_CUE,y_f_B_C./soilfactor,'LineWidth', linewidth)
ax = gca;
ax.XAxis.FontSize = fontsizes.xaxis;
ax.YAxis.FontSize = fontsizes.yaxis;
xlabel('days', 'FontSize', fontsizes.xlabel)
ylabel('amount carbon mg C g soil^{-1} per day ', 'FontSize', fontsizes.ylabel)

%set(gca,'ylim', [-2 1]);
legend(["R","R_O","f_U","f_BD_C", "f_Nec_C", "f_POM_C"], 'Interpreter', 'none', 'FontSize', fontsizes.legend) % , "f_B_growth"
%set(gca, 'FontSize', 14)
if title_on
    title(titel, 'Interpreter', 'none','FontSize', fontsizes.title)
end
set(gca,'TickDir',tick_dir);
set(figure10, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure10, fullfile(outputfolder, "fluxesC_" + appendix + ".png"));
savefig(fullfile(outputfolder, "fluxesC_" + appendix))





figure11 = figure('visible', isVisible) 
%plot(x_CUE,y_CUE,'LineWidth', linewidth)
%hold on 

plot(x_CUE,y_R_leaked./soilfactor,'LineWidth', linewidth, 'Color',colors.r)
hold on 
plot(x_CUE,y_f_N./soilfactor,'LineWidth', linewidth, 'Color',colors.dissolved_substrate)
hold on 
plot(x_CUE,y_f_BD_N./soilfactor,'LineWidth', linewidth, 'Color',colors.microbial_biomass)
hold on 
plot(x_CUE,y_f_MN_N./soilfactor,'LineWidth', linewidth, 'Color',colors.new_microbial_necromass)
hold on 
plot(x_CUE,y_f_POM_N./soilfactor,'LineWidth', linewidth, 'Color',colors.pom)
ax = gca;
ax.XAxis.FontSize = fontsizes.xaxis;
ax.YAxis.FontSize = fontsizes.yaxis;
xlabel('days', 'FontSize', fontsizes.xlabel)
ylabel('amount nitrogen mg N g soil^{-1} per day', 'FontSize', fontsizes.ylabel)

%set(gca,'ylim', [-2 1]);
legend(["R_leaked", "f_N", "f_BD_N", "f_Nec_N", "f_POM_N"], 'Interpreter', 'none', 'FontSize', fontsizes.legend)
%set(gca, 'FontSize', 14)
if title_on
    title(titel, 'Interpreter', 'none', 'FontSize', fontsizes.title)
end
set(gca,'TickDir',tick_dir);
set(figure11, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure11, fullfile(outputfolder, "fluxesN_" + appendix + ".png"));
savefig(fullfile(outputfolder, "fluxesN_" + appendix))




% figure12= figure('visible', isVisible) 
% growth = y_f_B_C - y_f_C + y_f_BD_C + y_R;
% phi = y_f_C./(y_C_S./y_N_S) - y_f_N;
% C_N_uptake =  y_f_C./y_f_N
% growth_N = -y_f_B_N + y_f_C/10 - y_f_BD_N;
% plot(x_CUE,phi,'LineWidth', linewidth)
% xlabel('days', 'FontSize', 14)
% ylabel('growth y_f_B_C ', 'FontSize', 14)
% %legend(["C_N"], 'Interpreter', 'none')
% set(gca, 'FontSize', 14)
% title(titel, 'Interpreter', 'none')
% set(figure12, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
% saveas(figure12, fullfile(outputfolder, "growth_" + appendix + ".png"));

figure13= figure('visible', isVisible) 
C_N = y_C_S./y_N_S;
plot(x_CUE,C_N,'LineWidth', linewidth, 'Color', colors.c_n)
ax = gca;
ax.XAxis.FontSize = fontsizes.xaxis;
ax.YAxis.FontSize = fontsizes.yaxis;
xlabel('days', 'FontSize', fontsizes.xlabel)
ylabel('C_N ', 'FontSize', fontsizes.ylabel)
legend(["C_N"], 'Interpreter', 'none', 'FontSize', fontsizes.legend)
%set(gca, 'FontSize', 14)
if title_on
    title(titel, 'Interpreter', 'none', 'FontSize', fontsizes.title)
end
set(gca,'TickDir',tick_dir);
set(figure13, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure13, fullfile(outputfolder, "C_N_" + appendix + ".png"));
savefig(fullfile(outputfolder, "C_N_" + appendix))
%pause

% figure14 = figure('visible', isVisible) 
% xlabel('days', 'FontSize', 14)
% plot(x_CUE,y_R_O,'LineWidth', linewidth)
% %plot(x_CUE,y_R_O./soilfactor,'LineWidth', linewidth)
% ylabel('amount nitrogen mg N g soil^{-1} per day ', 'FontSize', 14)
% 
% %set(gca,'ylim', [-2 1]);
% set(gca, 'FontSize', 14)
% title(titel, 'Interpreter', 'none')
% set(figure14, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
% saveas(figure14, fullfile(outputfolder, "phi_" + appendix + ".png"));

system("mogrify -trim " + outputfolder + "*.png")

close all











end


% figure
% C_N_B = y_C_B ./ y_N_B;
% C_N_S = y_C_S ./ y_N_S;
% C_N_MN = y_C_MN ./ y_N_MN;
% C_N_POM = y_C_POM ./ y_N_POM;
% plot(x, C_N_B, 'LineWidth', 2)
% hold on
% plot(x, C_N_S, 'LineWidth', 2)
% plot(x, C_N_MN, 'LineWidth', 2)
% plot(x, C_N_POM, 'LineWidth', 2)
% labels = {'C:N_B', 'C:N_S', 'C:N_MN', 'C:N_POM'};
% legend(labels, 'FontSize', 14)
% xlabel('days', 'FontSize', 14)
% ylabel('C:N ratio', 'FontSize', 14)
% set(gca, 'FontSize', 14)












% if false
% figure3 = figure('visible', isVisible) 
% linewidth = 4;
% plot(x,y_N_B./max(y_N_B),'LineWidth', linewidth)
% hold on
% plot(x,y_N_MN_old./max(y_N_MN_old),'LineWidth', linewidth)
% hold on
% plot(x,y_N_MN_new./max(y_N_MN_new),'LineWidth', linewidth)
% hold on
% plot(x,y_N_POM./max(y_N_POM),'LineWidth', linewidth)
% hold on
% plot(x,y_N_S./max(y_N_S),'LineWidth', linewidth)
% hold on
% plot(x,y_sumleakedC_N./max(y_sumleakedC_N),'LineWidth', linewidth)
% legend(names, 'FontSize', 14, 'Location','southeast')
% xlabel('days', 'FontSize', 14)
% ylabel('amount nitrogen normalized by the maximum value', 'FontSize', 14)
% ylim([-0.1,1.1])
% set(gca, 'FontSize', 14)
% title(titel, 'Interpreter', 'none')
% set(figure3, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
% saveas(figure3, fullfile(outputfolder, "C_norm_" + appendix + ".png"));
% 
% figure4 = figure('visible', isVisible) 
% linewidth = 4;
% plot(x,y_C_B./max(y_C_B),'LineWidth', linewidth)
% hold on
% plot(x,y_C_MN_old./max(y_C_MN_old),'LineWidth', linewidth)
% hold on
% plot(x,y_C_MN_new./max(y_C_MN_new),'LineWidth', linewidth)
% hold on
% plot(x,y_C_POM./max(y_C_POM),'LineWidth', linewidth)
% hold on
% plot(x,y_C_S./max(y_C_S),'LineWidth', linewidth)
% hold on
% plot(x,y_CO2./max(y_CO2),'LineWidth', linewidth)
% legend(names, 'FontSize', 14, 'Location','southeast')
% xlabel('days', 'FontSize', 14)
% ylabel('amount carbon normalized by the maximum value', 'FontSize', 14)
% ylim([-0.1,1.1])
% set(gca, 'FontSize', 14)
% title(titel, 'Interpreter', 'none')
% set(figure4, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
% saveas(figure4, fullfile(outputfolder, "N_norm_last_" + appendix + ".png"));
% %print(fullfile(outputfolder, 'figure_large_dpi.png'), '-dpng', '-r300');
% 
% 
% 
% 
% 
% figure6 = figure('visible', isVisible) 
% linewidth = 4;
% plot(x,y_N_B./max(y_N_B),'LineWidth', linewidth)
% hold on
% plot(x,y_N_MN_old./y_N_MN_old(end),'LineWidth', linewidth)
% hold on
% plot(x,y_N_MN_new./y_N_MN_new(end),'LineWidth', linewidth)
% hold on
% plot(x,y_N_POM./y_N_POM(end),'LineWidth', linewidth)
% hold on
% plot(x,y_N_S./max(y_N_S),'LineWidth', linewidth)
% hold on
% plot(x,y_sumleakedC_N./y_sumleakedC_N(end),'LineWidth', linewidth)
% legend(names, 'FontSize', 14, 'Location','southeast')
% xlabel('days', 'FontSize', 14)
% ylabel('amount nitrogen normalized by the steady state value', 'FontSize', 14)
% %ylim([-0.1,1.1])
% set(gca, 'FontSize', 14)
% title(titel, 'Interpreter', 'none')
% set(figure6, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
% saveas(figure6, fullfile(outputfolder, "C_norm_last_" + appendix + ".png"));
% 
% figure7 = figure('visible', isVisible) 
% linewidth = 4;
% plot(x,y_C_B./max(y_C_B),'LineWidth', linewidth)
% hold on
% plot(x,y_C_MN_old./y_C_MN_old(end),'LineWidth', linewidth)
% hold on
% plot(x,y_C_MN_new./y_C_MN_new(end),'LineWidth', linewidth)
% hold on
% plot(x,y_C_POM./y_C_POM(end),'LineWidth', linewidth)
% hold on
% plot(x,y_C_S./max(y_C_S),'LineWidth', linewidth)
% hold on
% plot(x,y_CO2./y_CO2(end),'LineWidth', linewidth)
% legend(names, 'FontSize', 14, 'Location','southeast')
% xlabel('days', 'FontSize', 14)
% ylabel('amount carbon normalized by the steady state value', 'FontSize', 14)
% %ylim([-0.1,1.1])
% set(gca, 'FontSize', 14)
% title(titel, 'Interpreter', 'none')
% set(figure7, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
% saveas(figure7, fullfile(outputfolder, "N_norm_" + appendix + ".png"));
% %print(fullfile(outputfolder, 'figure_large_dpi.png'), '-dpng', '-r300');
% end