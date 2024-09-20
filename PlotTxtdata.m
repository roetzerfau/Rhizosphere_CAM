function PlotTxtdata()
close all
figure
root = "/home.local/roetzer/C_N/txtdata/";
%root = "txtdata_move/";% (copy)/";_mass_balance2
percentage_plot = 25;


fileID = fopen(root + 'C_BVector.txt','r');
formatSpec = '%f %f';
C_B = fscanf(fileID,formatSpec,  [2 Inf])';
y_C_B = C_B(:,2);
x_C_B = C_B(:,1);

%plot(x_C_B, y_C_B, 'DisplayName', 'C_B');
%xlabel('days')
%ylabel('carbon biomass')
%hold on 

fileID = fopen(root + 'C_SVector.txt','r');
formatSpec = '%f %f';
C_S = fscanf(fileID,formatSpec,  [2 Inf])';
y_C_S = C_S(:,2);
x_C_S = C_S(:,1);
%plot(x_C_S, y_C_S, 'DisplayName', 'C_S');

%hold on 

fileID = fopen( root + 'C_MNVector.txt','r');
formatSpec = '%f %f';
C_MN = fscanf(fileID,formatSpec,  [2 Inf])';
y_C_MN = C_MN(:,2);
x_C_MN = C_MN(:,1);
%plot(x_MN, y_MN, 'DisplayName', 'MN');
%xlabel('days')
%ylabel('amount necromass')

fileID = fopen( root + 'C_POMconcVector.txt','r');
formatSpec = '%f %f';
C_POM = fscanf(fileID,formatSpec,  [2 Inf])';
y_C_POM = C_POM(:,2);
x_C_POM = C_POM(:,1);

fileID = fopen( root + 'CO2Vector.txt','r');
formatSpec = '%f %f';
CO2 = fscanf(fileID,formatSpec,  [2 Inf])';
y_CO2 = CO2(:,2);
x_CO2 = CO2(:,1);

x = x_C_MN;
y = [];
for i = 1:numel(y_C_MN)
   y_add = [y_C_B(i), y_C_S(i),y_C_MN(i), y_C_POM(i), y_CO2(i)];
   sum_C = y_C_B(i) + y_C_S(i) + y_C_MN(i) + y_C_POM(i) + y_CO2(i)
   y_add = y_add/(250*250);
   y = [y; y_add];
end
y
names = {'biomass', 'dissolved substrate', 'Necromass', 'POM', 'CO2'};
bar(x(1:end),y,'stacked')
legend(names)
xlabel('days')
ylabel('amount carbon g cm^-3')

% Add percentage annotations
for k = 1:size(y, 1)
    for j = 1:size(y, 2)
        if(mod(k,percentage_plot)== 0)
            percentage = y(k, j) / sum(y(k, :)) * 100;
            text(x(k), sum(y(k, 1:j)) - y(k, j)/2, sprintf('%i%', ceil(percentage)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle')
        end
    end
end
%.1f






figure






fileID = fopen(root + 'N_BVector.txt','r');
formatSpec = '%f %f';
N_B = fscanf(fileID,formatSpec,  [2 Inf])';
y_N_B = N_B(:,2);
x_N_B = N_B(:,1);

%plot(x_C_B, y_C_B, 'DisplayName', 'C_B');
%xlabel('days')
%ylabel('carbon biomass')
%hold on 

fileID = fopen(root + 'N_SVector.txt','r');
formatSpec = '%f %f';
N_S = fscanf(fileID,formatSpec,  [2 Inf])';
y_N_S = N_S(:,2);
x_N_S = N_S(:,1);
%plot(x_C_S, y_C_S, 'DisplayName', 'C_S');

%hold on 

fileID = fopen( root + 'N_MNVector.txt','r');
formatSpec = '%f %f';
C_MN = fscanf(fileID,formatSpec,  [2 Inf])';
y_N_MN = C_MN(:,2);
x_N_MN = C_MN(:,1);
%plot(x_MN, y_MN, 'DisplayName', 'MN');
%xlabel('days')
%ylabel('amount necromass')

fileID = fopen( root + 'N_POMconcVector.txt','r');
formatSpec = '%f %f';
N_POM = fscanf(fileID,formatSpec,  [2 Inf])';
y_N_POM = N_POM(:,2);
x_N_POM = N_POM(:,1);

fileID = fopen( root + 'sumleakedN_S.txt','r');
formatSpec = '%f %f';
sumleakedC_N = fscanf(fileID,formatSpec,  [2 Inf])';
y_sumleakedC_N = sumleakedC_N(:,2);
x_sumleakedC_N = sumleakedC_N(:,1);

x = x_N_MN;
y = [];

for i = 1:numel(y_N_MN)
   y_add = [y_N_B(i), y_N_S(i),y_N_MN(i), y_N_POM(i), y_sumleakedC_N(i)];
   sum_N = y_N_B(i) + y_N_S(i) + y_N_MN(i) + y_N_POM(i) + y_sumleakedC_N(i)
   y_add = y_add/(250*250);
   y = [y; y_add];
end
y
names = {'biomass', 'dissolved substrate', 'Necromass', 'POM', 'leaked N'};
bar(x(1:end),y,'stacked')
legend(names)
xlabel('days')
ylabel('amount nitrogen g cm^-3 ')
% Add percentage annotations
for k = 1:size(y, 1)
    for j = 1:size(y, 2)
        if(mod(k,percentage_plot)== 0)
            percentage = y(k, j) / sum(y(k, :)) * 100;
            text(x(k), sum(y(k, 1:j)) - y(k, j)/2, sprintf('%i%', ceil(percentage)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle')
        end
    end
end
%.1f


figure
C_N_B = y_C_B./y_N_B;
C_N_S = y_C_S./y_N_S;
C_N_MN = y_C_MN./y_N_MN;
C_N_POM = y_C_POM./y_N_POM;
plot(x,C_N_B)
hold on
plot(x,C_N_S)
hold on
plot(x,C_N_MN)
hold on
plot(x,C_N_POM)

labels = {'C:N_B', 'C:N_S','C:N_MN','C:N_POM' };
legend(labels)

% figure 
% for i = 1:numel(y_N_MN)
%    sum_C(i) = y_C_B(i) + y_C_S(i) + y_C_MN(i) + y_C_POM(i);
%    sum_N(i) = y_N_B(i) + y_N_S(i) + y_N_MN(i) + y_N_POM(i);
% end
% C_N_total = sum_C./sum_N;
% plot(x,C_N_total)
% labels = {'C:N_total' };
% legend(labels)
end