clear; close all;

figure
step = 30;
ende = 330;
%% agg33
a = load('/home.other/fauam9/roetzer/Frontier22_maxDecay/agg34_noMuc/Rhizosphere_CAM/FinalConfig/config.5.mat');
b = load('/home.other/fauam9/roetzer/Frontier22_maxDecay/agg34_noMuc/Rhizosphere_CAM/FinalConfig/config.100.mat');

b_root = load('/home.other/fauam9/roetzer/Frontier22_maxDecay/agg34_noMuc/Rhizosphere_CAM/FinalConfig/rootConfig.100.mat');
a_root = load('/home.other/fauam9/roetzer/Frontier22_maxDecay/agg34_noMuc/Rhizosphere_CAM/FinalConfig/rootConfig.5.mat');


r_a = sqrt(sum(a_root.rootVector)/pi)*2;
r_b = sqrt(sum(b_root.rootVector)/pi)*2;
r = abs(r_b - r_a);

disp = calculateDisplacementVectors(a.g, a_root.rootVector, b.particleList, a.particleList);
U = disp(:,1)*2;
V = disp(:,2)*2;

dist = disp(:,5) * 2;

M = sqrt(U.^2 + V.^2);

d_bar = M/r * 100;
edges = 0:step:ende;
%d_bar_dis = edges;
Y = discretize(dist,edges);
edges = edges(2:end);
for index = 1:(numel(edges))
    d_bar_dis(index) = mean(d_bar(Y == index));
    %min = min(d_bar(Y == index));
    S(index) = std(M(Y == index)/r *100);
end
x1 = (edges );%+ step/2
y1 = d_bar_dis;
errlow1 =  S;
errhigh1 =  S;

%% agg 18

a = load('/home.other/fauam9/roetzer/Frontier22_maxDecay/agg19_noMuc/Rhizosphere_CAM/FinalConfig/config.5.mat');
b = load('/home.other/fauam9/roetzer/Frontier22_maxDecay/agg19_noMuc/Rhizosphere_CAM/FinalConfig/config.100.mat');

b_root = load('/home.other/fauam9/roetzer/Frontier22_maxDecay/agg19_noMuc/Rhizosphere_CAM/FinalConfig/rootConfig.100.mat');
a_root = load('/home.other/fauam9/roetzer/Frontier22_maxDecay/agg19_noMuc/Rhizosphere_CAM/FinalConfig/rootConfig.5.mat');


r_a = sqrt(sum(a_root.rootVector)/pi)*2;
r_b = sqrt(sum(b_root.rootVector)/pi)*2;
r = abs(r_b - r_a);

disp = calculateDisplacementVectors(a.g, a_root.rootVector, b.particleList, a.particleList);
U = disp(:,1)*2;
V = disp(:,2)*2;

dist = disp(:,5)*2;

M = sqrt(U.^2 + V.^2);

d_bar = M/r * 100;
edges = 0:step:ende;

%d_bar_dis = edges;
Y = discretize(dist,edges);
edges = edges(2:end);
for index = 1:(numel(edges))
    d_bar_dis(index) = mean(d_bar(Y == index));

    S(index) = std(M(Y == index)/r *100);
end
x2 = (edges );
y2 = d_bar_dis;
errlow2 =  S;
errhigh2 =  S;







%% display
%colormap jet
y = [y1;y2]';
err = [errlow1;errlow2]';
b = bar([x1; x2]',[y1;y2]','grouped','BarWidth',1 );
b(2).FaceColor = [0.8500 0.3250 0.0980];
b(1).FaceColor = [0.9290 0.6940 0.1250];




% hold on
% ngroups = size(y, 1);
% nbars = size(y, 2);
% %Calculating the width for each bar group
% groupwidth =10;% min(0.8, nbars/(nbars + 1.5 ));
% for i = 1:nbars
%     x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
%     errorbar(x, y(:,i), err(:,i), '.');
% end
% hold off



% hold on
% errorbar([x1; x2]',[y1;y2]',[errlow1;errlow2]',[errhigh1;errhigh2]');
% er.Color = [0 0 0];                            
% er.LineStyle = 'none';
% hold off


size = 20;
set(gca,"FontSize",size-5)
%title('Displacement of solid building units around growing root')
c ="Distance from root ("+char(181) + "m)";
xlabel(c, 'FontSize', size)
ylabel('Displacement normed by root expansion (%)', 'FontSize', size)
xticks(0:25*2:ende);
label = ["soil with 33% clay", "soil with 18% clay"] 
legend(label, 'FontSize', size)
%imwrite(getframe(gcf).cdata, 'displacementEval.png')
saveas(gcf, 'displacementEval.png')
%figure 
%geoBulk = reshape(b.bulkVector, [N,N]);
%imshow(geoBulk)
