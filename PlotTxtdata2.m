function PlotTxtdata2(file, description)
close all
if nargin < 2
    description =file;
end
if nargin < 1
    file = "noMove_longterm";
end
imageFolder = "bilder/CN/";
appendix = file%"noMove_longterm_noreactiveEdgesforDecay" %"noMove_longterm" %mucilageC_after100stepsCN01
outputfolder = imageFolder + appendix +"/";
if ~exist(outputfolder, 'dir')  % Check if the folder exists
    mkdir(outputfolder);        % Create the folder
    disp(['Folder created at: ', outputfolder]);
else
    disp(['Folder already exists at: ', outputfolder]);
end
%ANNZ / SNNZ    * (sum(C_A)/ANNZ)*  (C_SNNZ/sum(C_S)  
%fläche anteil/ fläche soil   *  anteil densitiy/solid density 
porosity = 0.45;
soil_particleNNZ = 250 * 250 * (1-porosity);
soilParticleDensity = 2.65;
isNNZ = false;

root = "/home.local/roetzer/C_N/tx2tdata_" + appendix +"/";
titel = description;
percentage_plot = 50;%;%25;
isVisible = 'on';
nof_cells = 1;
if file.contains("2years")
    last_value = 2000;
else
last_value = Inf;
end
% DOC ist zu hoch

ylim_C_stacked = [0 0.2];
ylim_N_stacked = [0 4*10^-3];

ylim_C = [10^-8 10^-0];
ylim_N = [10^-8 10^0 ];
%% carbon data
fileID = fopen(root + 'C_BVector.txt','r');
formatSpec = '%f %f';
C_B = fscanf(fileID,formatSpec,  [2 last_value])';
y_C_B = C_B(:,2)/nof_cells;
x_C_B = C_B(:,1);

if(isNNZ)
fileID = fopen(root + 'N_BVectorNNZ.txt','r');
formatSpec = '%f %f';
C_B = fscanf(fileID,formatSpec,  [2 last_value])';
y_C_BNNZ = C_B(:,2);
x_C_BNNZ = C_B(:,1);
end

fileID = fopen(root + 'C_SVector.txt','r');
formatSpec = '%f %f';
C_S = fscanf(fileID,formatSpec,  [2 last_value])';
y_C_S = C_S(:,2)/nof_cells;
x_C_S = C_S(:,1);

if(isNNZ)
fileID = fopen(root + 'C_SVectorNNZ.txt','r');
formatSpec = '%f %f';
C_S = fscanf(fileID,formatSpec,  [2 last_value])';
y_C_SNNZ = C_S(:,2);
x_C_SNNZ = C_S(:,1);
end

fileID = fopen( root + 'C_MNVector.txt','r');
formatSpec = '%f %f';
C_MN = fscanf(fileID,formatSpec,  [2 last_value])';
y_C_MN = C_MN(:,2)/nof_cells;
x_C_MN = C_MN(:,1);


if(isNNZ)
fileID = fopen( root + 'C_MNVectorNNZ.txt','r');
formatSpec = '%f %f';
C_MN = fscanf(fileID,formatSpec,  [2 last_value])';
y_C_MNNNZ = C_MN(:,2);
x_C_MNNNZ = C_MN(:,1);
end

fileID = fopen( root + 'C_MNVector_old.txt','r');
formatSpec = '%f %f';
C_MN_old = fscanf(fileID,formatSpec,  [2 last_value])';
y_C_MN_old = C_MN_old(:,2)/nof_cells;
x_C_MN_old = C_MN_old(:,1);

if(isNNZ)
fileID = fopen( root + 'C_MNVector_oldNNZ.txt','r');
formatSpec = '%f %f';
C_MN_old = fscanf(fileID,formatSpec,  [2 last_value])';
y_C_MN_oldNNZ = C_MN_old(:,2);
x_C_MN_oldNNZ = C_MN_old(:,1);
end

fileID = fopen( root + 'C_MNVector_new.txt','r');
formatSpec = '%f %f';
C_MN_new = fscanf(fileID,formatSpec,  [2 last_value])';
y_C_MN_new = C_MN_new(:,2)/nof_cells;
x_C_MN_new = C_MN_new(:,1);

if(isNNZ)
fileID = fopen( root + 'C_MNVector_newNNZ.txt','r');
formatSpec = '%f %f';
C_MN_new = fscanf(fileID,formatSpec,  [2 last_value])';
y_C_MN_newNNZ = C_MN_new(:,2);
x_C_MN_newNNZ = C_MN_new(:,1);
end

fileID = fopen( root + 'C_POMconcVector.txt','r');
formatSpec = '%f %f';
C_POM = fscanf(fileID,formatSpec,  [2 last_value])';
y_C_POM = C_POM(:,2)/nof_cells;
x_C_POM = C_POM(:,1);

if(isNNZ)
fileID = fopen( root + 'C_POMconcVectorNNZ.txt','r');
formatSpec = '%f %f';
C_POM = fscanf(fileID,formatSpec, [2 last_value])';
y_C_POMNNZ = C_POM(:,2);
x_C_POMNNZ = C_POM(:,1);
end

fileID = fopen( root + 'CO2Vector.txt','r');
formatSpec = '%f %f';
CO2 = fscanf(fileID,formatSpec,  [2 last_value])';
y_CO2 = CO2(:,2)/nof_cells;
x_CO2 = CO2(:,1);

if(isNNZ)
fileID = fopen( root + 'CO2VectorNNZ.txt','r');
formatSpec = '%f %f';
CO2 = fscanf(fileID,formatSpec,  [2 last_value])';
y_CO2NNZ = CO2(:,2);
x_CO2NNZ = CO2(:,1);
end

x = x_C_MN;
y_C = [];
y_C_sum = [];
CO2_days = [];
for i = 1:numel(y_C_MN)
   y_C_add = [y_C_B(i), y_C_MN_old(i),y_C_MN_new(i),y_C_POM(i),  y_C_S(i),y_CO2(i)];
   %y_C_add = [y_C_B(i)/y_C_BNNZ(i), y_C_S(i)/y_C_SNNZ(i),y_C_MN_old(i)/y_C_MN_oldNNZ(i),y_C_MN_new(i)/y_C_MN_newNNZ(i), y_C_POM(i)/y_C_POMNNZ(i), y_CO2(i)/y_CO2NNZ(i)];
   sum_C = y_C_B(i)+ y_C_MN_old(i)+y_C_MN_new(i)+y_C_POM(i)+  y_C_S(i);
   y_C_add = y_C_add./( soil_particleNNZ * soilParticleDensity/1000 );% / soilParticleDensity;
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
TOC = y_C_sum./(soil_particleNNZ * soilParticleDensity/1000); 
CO2_hour_avag = mean(CO2_days);

%% Nitrogen data
%indices = find(POMageVector > 0 & POMageVector ~= max(POMageVector))
fileID = fopen(root + 'N_BVector.txt','r');
formatSpec = '%f %f';
N_B = fscanf(fileID,formatSpec,  [2 last_value])';
y_N_B = N_B(:,2)/nof_cells;
x_N_B = N_B(:,1);

if(isNNZ)
fileID = fopen(root + 'N_BVectorNNZ.txt','r');
formatSpec = '%f %f';
N_B = fscanf(fileID,formatSpec,  [2 last_value])';
y_N_BNNZ = N_B(:,2);
x_N_BNNZ = N_B(:,1);
end

fileID = fopen(root + 'N_SVector.txt','r');%NNZ
formatSpec = '%f %f';
N_S = fscanf(fileID,formatSpec,  [2 last_value])';
y_N_S = N_S(:,2)/nof_cells;
x_N_S = N_S(:,1);

if(isNNZ)
fileID = fopen(root + 'N_SVectorNNZ.txt','r');
formatSpec = '%f %f';
N_S = fscanf(fileID,formatSpec,  [2 last_value])';
y_N_SNNZ = N_S(:,2);
x_N_SNNZ = N_S(:,1);
end

fileID = fopen( root + 'N_MNVector.txt','r');
formatSpec = '%f %f';
N_MN = fscanf(fileID,formatSpec,  [2 last_value])';
y_N_MN = N_MN(:,2)/nof_cells;
x_N_MN = N_MN(:,1);

if(isNNZ)
fileID = fopen( root + 'N_MNVectorNNZ.txt','r');
formatSpec = '%f %f';
N_MNNNZ = fscanf(fileID,formatSpec,  [2 last_value])';
y_N_MNNNZ = N_MNNNZ(:,2);
x_N_MNNNZ = N_MNNNZ(:,1);
end

fileID = fopen( root + 'N_MNVector_old.txt','r');
formatSpec = '%f %f';
N_MN_old = fscanf(fileID,formatSpec,  [2 last_value])';
y_N_MN_old = N_MN_old(:,2)/nof_cells;
x_N_MN_old = N_MN_old(:,1);

if(isNNZ)
fileID = fopen( root + 'N_MNVector_oldNNZ.txt','r');
formatSpec = '%f %f';
N_MN_old = fscanf(fileID,formatSpec,  [2 last_value])';
y_N_MN_oldNNZ = N_MN_old(:,2);
x_N_MN_oldNNZ = N_MN_old(:,1);
end

fileID = fopen( root + 'N_MNVector_new.txt','r');
formatSpec = '%f %f';
N_MN_new = fscanf(fileID,formatSpec,  [2 last_value])';
y_N_MN_new = N_MN_new(:,2)/nof_cells;
x_N_MN_new = N_MN_new(:,1);

if(isNNZ)
fileID = fopen( root + 'N_MNVector_newNNZ.txt','r');
formatSpec = '%f %f';
N_MN_new = fscanf(fileID,formatSpec,  [2 last_value])';
y_N_MN_newNNZ = N_MN_new(:,2);
x_N_MN_newNNZ = N_MN_new(:,1);
end

fileID = fopen( root + 'N_POMconcVector.txt','r');
formatSpec = '%f %f';
N_POM = fscanf(fileID,formatSpec,  [2 last_value])';
y_N_POM = N_POM(:,2)/nof_cells;
x_N_POM = N_POM(:,1);

if(isNNZ)
fileID = fopen( root + 'N_POMconcVectorNNZ.txt','r');
formatSpec = '%f %f';
N_POM = fscanf(fileID,formatSpec,  [2 last_value])';
y_N_POMNNZ = N_POM(:,2);
x_N_POMNNZ = N_POM(:,1);
end


fileID = fopen( root + 'leakedNVector.txt','r');
formatSpec = '%f %f';
leakedN = fscanf(fileID,formatSpec,  [2 Inf])';
y_sumleakedC_N = leakedN(:,2)/nof_cells;
x_leakedN = leakedN(:,1);

if(isNNZ)
fileID = fopen( root + 'leakedNVectorNNZ.txt','r');
formatSpec = '%f %f';
leakedNNNZ = fscanf(fileID,formatSpec,  [2 Inf])';
y_leakedNNNZ = leakedNNNZ(:,2);
x_leakedNNNZ = leakedNNNZ(:,1);
end

% fileID = fopen( root + 'sumleakedN_S.txt','r');
% formatSpec = '%f %f';
% sumleakedC_N = fscanf(fileID,formatSpec,  [2 last_value])';
% y_sumleakedC_N = sumleakedC_N(:,2)/nof_cells;
% x_sumleakedC_N = sumleakedC_N(:,1);

x = x_N_MN;
y_N = [];
for i = 1:numel(y_N_MN)
   %y_add = [y_N_B(i), y_N_S(i), y_N_MN(i), y_N_POM(i), y_sumleakedC_N(i)];
   y_N_add = [y_N_B(i), y_N_MN_old(i),y_N_MN_new(i), y_N_POM(i),y_N_S(i),  y_sumleakedC_N(i)];
   %y_N_add = [y_N_B(i)/y_N_BNNZ(i), y_N_MN_old(i)/y_N_MN_oldNNZ(i),y_N_MN_new(i)/y_N_MN_newNNZ(i), y_N_POM(i)/y_N_POMNNZ(i),y_N_S(i)/y_N_SNNZ(i),  y_leakedN(i)/y_leakedNNNZ(i)];
   %sum_N = y_N_B(i) + y_N_S(i) + y_N_MN(i) + y_N_POM(i) + y_sumleakedC_N(i);
   y_N_add = y_N_add./( soil_particleNNZ * soilParticleDensity/1000 );%/ soilParticleDensity;
   y_N = [y_N; y_N_add];
end

%% Plots 
percentage_plot = round(numel(x_N_POM)/20);%;%25;
linewidth = 4;
if false

%% figure 1
figure1 = figure('visible', isVisible)
names = {'biomass', 'dissolved substrate', 'Necromass', 'POM', 'CO2'};
names = {'biomass',  'Necromass old', 'Necromass new', 'POM', 'dissolved substrate','CO2'};
bar(x(1:end), y_C, 'stacked')
legend(names, 'FontSize', 14, 'Location','southeast')
xlabel('days', 'FontSize', 14)
%ylabel('amount carbon g cm^{-3}', 'FontSize', 14)
ylabel('amount carbon mg C g soil^{-1}', 'FontSize', 14)
%set(gca,'ylim',ylim_C_stacked);
title(titel , 'Interpreter', 'none')
set(gca, 'FontSize', 14)

% Add percentage annotations
for k = 1:size(y_C, 1)
    for j = 1:size(y_C, 2)
        if mod(k, percentage_plot) == 0
            percentage = y_C(k, j) / sum(y_C(k, :)) * 100;
            text(x(k), sum(y_C(k, 1:j)) - y_C(k, j)/2, sprintf('%i%%', ceil(percentage)), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 12)
        end
    end
end
set(figure1, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure1, fullfile(outputfolder, "C_dist_" + appendix + ".png"));




%% figure 2
figure2 = figure('visible', isVisible)
names = {'biomass', 'dissolved substrate', 'Necromass', 'POM', 'leaked N'};
names = {'biomass', 'Necromass old','Necromass new', 'POM', 'dissolved substrate', 'leaked N'};
bar(x(1:end), y_N, 'stacked')
legend(names, 'FontSize', 14, 'Location','southeast')
xlabel('days', 'FontSize', 14)
%ylabel('amount nitrogen g cm^{-3}', 'FontSize', 14)
ylabel('amount nitrogen mg N g soil^{-1}', 'FontSize', 14)
title(titel, 'Interpreter', 'none')
%set(gca,'ylim',ylim_N_stacked);
set(gca, 'FontSize', 14)

% Add percentage annotations
for k = 1:size(y_N, 1)
    for j = 1:size(y_N, 2)
        if mod(k, percentage_plot) == 0
            percentage = y_N(k, j) / sum(y_N(k, :)) * 100;
            text(x(k), sum(y_N(k, 1:j)) - y_N(k, j)/2, sprintf('%i%%', ceil(percentage)), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 12)
        end
    end
end
set(figure2, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure2, fullfile(outputfolder, "N_dist_" + appendix + ".png"));

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
names = {'biomass', 'Necromass old','Necromass new', 'POM', 'dissolved substrate', 'leaked N'};
linewidth = 4;
plot(x,y_N_B,'LineWidth', linewidth)
hold on
plot(x,y_N_MN_old,'LineWidth', linewidth)
hold on
plot(x,y_N_MN_new,'LineWidth', linewidth)
hold on
plot(x,y_N_POM,'LineWidth', linewidth)
hold on
plot(x,y_N_S,'LineWidth', linewidth)
hold on
plot(x,y_sumleakedC_N,'LineWidth', linewidth)
legend(names, 'FontSize', 14, 'Location','southeast')
xlabel('days', 'FontSize', 14)
%ylabel('amount nitrogen g cm^{-3} (log scale)', 'FontSize', 14)
ylabel('amount nitrogen mg N g soil^{-1} (log scale)', 'FontSize', 14)
% if(~(file.contains("normal") || file.contains("noPOMDecay")))
% set(gca,'ylim',ylim_N);
% end
set(gca, 'YScale', 'log');  % Set y-axis to log scale


%ylim = mylimits;
set(gca, 'FontSize', 14)
title(titel, 'Interpreter', 'none')
set(figure5, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure5, fullfile(outputfolder, "N_log_" + appendix + ".png"));

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


plot(x,y_C_B,'LineWidth', linewidth)
hold on
plot(x,y_C_MN_old,'LineWidth', linewidth)
hold on
plot(x,y_C_MN_new,'LineWidth', linewidth)
hold on
plot(x,y_C_POM,'LineWidth', linewidth)
hold on
plot(x,y_C_S,'LineWidth', linewidth)
hold on
plot(x,y_CO2,'LineWidth', linewidth)

names = {'biomass',  'Necromass old', 'Necromass new', 'POM', 'dissolved substrate','CO2'};
legend(names, 'FontSize', 14, 'Location','southeast')
xlabel('days', 'FontSize', 14)
%ylim(mylimits);
%ylabel('amount carbon g cm^{-3} (log scale)', 'FontSize', 14)
ylabel('amount carbon mg C g soil^{-1} (log scale)', 'FontSize', 14)
% if(~(file.contains("normal") || file.contains("noPOMDecay")))
% set(gca,'ylim',ylim_C);
% end
set(gca, 'YScale', 'log');  % Set y-axis to log scale
%ylim([-0.1,1.1])
set(gca, 'FontSize', 14)
title(titel, 'Interpreter', 'none')
set(figure6, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure6, fullfile(outputfolder, "C_log_" + appendix + ".png"));
%print(fullfile(outputfolder, 'figure_large_dpi.png'), '-dpng', '-r300');




end

fileID = fopen( root + 'CUE.txt','r');
formatSpec = '%f %f';
CUE = fscanf(fileID,formatSpec,  [2 last_value])';
y_CUE = CUE(:,2);
x_CUE = CUE(:,1);


fileID = fopen( root + 'B.txt','r');
formatSpec = '%f %f';
B = fscanf(fileID,formatSpec,  [2 last_value])';
y_B = B(:,2);
x_B = B(:,1);

fileID = fopen( root + 'R.txt','r');
formatSpec = '%f %f';
R = fscanf(fileID,formatSpec,  [2 last_value])';
y_R = R(:,2);
x_R = R(:,1);

fileID = fopen( root + 'T.txt','r');
formatSpec = '%f %f';
T = fscanf(fileID,formatSpec,  [2 last_value])';
y_T = T(:,2);
x_T = T(:,1);

fileID = fopen( root + 'U.txt','r');
formatSpec = '%f %f';
U = fscanf(fileID,formatSpec,  [2 last_value])';
y_U = U(:,2);
x_U = U(:,1);










figure8 = figure('visible', isVisible) 
%plot(x_CUE,y_CUE,'LineWidth', linewidth)
%hold on 
CUE_A = 1 - y_R./y_U;
plot(x_CUE,CUE_A,'LineWidth', linewidth)
hold on 
CUE_B = 1 - (y_R+ y_T)./y_U;
plot(x_CUE,CUE_B,'LineWidth', linewidth)
xlabel('days', 'FontSize', 14)
ylabel('Carbon Use Efficiency Microbes', 'FontSize', 14)
set(gca,'ylim', [-2 1]);
legend(["CUE_A", "CUE_B"])
set(gca, 'FontSize', 14)
title(titel, 'Interpreter', 'none')
set(figure8, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure8, fullfile(outputfolder, "CUE_" + appendix + ".png"));












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