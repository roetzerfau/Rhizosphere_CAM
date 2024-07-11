function PlotTxtdata()
close all
figure

fileID = fopen('/home.local/roetzer/09_06/txtdata/C_BVector.txt','r');
formatSpec = '%f %f';
C_B = fscanf(fileID,formatSpec,  [2 Inf])';
y_C_B = C_B(:,2);
x_C_B = C_B(:,1);

plot(x_C_B, y_C_B, 'DisplayName', 'C_B');
xlabel('days')
ylabel('carbon biomass')
%hold on 

% fileID = fopen('txtdata_move/C_SVector.txt','r');
% formatSpec = '%f %f';
% C_S = fscanf(fileID,formatSpec,  [2 Inf])';
% y_C_S = C_S(:,2)
% x_C_S = C_S(:,1)
% 
% 
% plot(x_C_S, y_C_S, 'DisplayName', 'C_S');
% hold on 

fileID = fopen('/home.local/roetzer/09_06/txtdata/MNVector.txt','r');
formatSpec = '%f %f';
MN = fscanf(fileID,formatSpec,  [2 Inf])';
y_MN = MN(:,2);
x_MN = MN(:,1);
plot(x_MN, y_MN, 'DisplayName', 'MN');
xlabel('days')
ylabel('amount necromass')

end