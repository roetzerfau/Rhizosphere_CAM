clear; close all;
a = load('FinalConfig/config.5.mat');
b = load('FinalConfig/config.100.mat');
N = a.g.NX;

b_root = load('FinalConfig/rootConfig.100.mat');
b_root_geo = reshape(b_root.rootVector, [N,N]);
%figure
%imshow(b_root_geo)
a_root = load('FinalConfig/rootConfig.5.mat');

c =  b.particleList(1:2);

r_a = sqrt(sum(a_root.rootVector)/pi);
r_b = sqrt(sum(b_root.rootVector)/pi);
r = abs(r_b - r_a);

disp = calculateDisplacementVectors(a.g, a_root.rootVector, b.particleList(1:end), a.particleList(1:end));

geoBulk_color = a.bulkVector;
geoBulk = reshape(a.bulkVector, [N,N]);
%geoBulk = cat(3, geoBulk, geoBulk, geoBulk);


 
X = disp(:,3);
Y = disp(:,4);
U = disp(:,1);
V = disp(:,2);
%V(:) = 0;
dist = disp(:,5);
ind = sub2ind(size(geoBulk), X, Y);
%ind = a.particleList{1};
[row,col] = ind2sub(size(geoBulk),ind);
% row = X(1);
% col = Y(2);
%geoBulk(row,col) = 0.5;
% geoBulk(row,col,1) = 1;
% geoBulk(row,col,2:3) = 0;
geoBulk(ind) = 0.5;

%figure
%imshow(~geoBulk, [[0.75,0.75,0.75];[1,1,1]])
%hold on
%quiver(Y,X,V,U, 0, 'LineWidth',2)

%imwrite(getframe(gcf).cdata, 'displacementParticle.png')%,"ScreenSize",[1000 1060]
%saveas(double(gcf), 'displacementParticle.png')
M = sqrt(U.^2 + V.^2);

figure
for i = 1:numel(M)
particleInd = a.particleList{i};
geoBulk_color(particleInd) = geoBulk_color(particleInd) .* (M(i) + 1);
end
max_value = max(geoBulk_color);
%geoBulk_color = geoBulk_color ./(max_value * 255);
geoBulk_color = reshape(geoBulk_color, [N,N]);
geoBulk_color = geoBulk_color(1:end,50:150);
imshow(geoBulk_color, [])

%map = cool;
%map(1,1:3) = 1;
map = flipud(hot);
colormap(map)
colorbar

d_bar = M/r * 100;
edges = 0:25:200;
d_bar_dis = edges;
Y = discretize(dist,edges);
for index = 1:numel(edges)
    d_bar_dis(index) = mean(d_bar(Y == index));
end

% [B,I] = sort(dist);
% M = d_bar(I);
figure
%scatter(dist, d_bar);
bar(edges + 12.5,d_bar_dis)
title('Displacement of solid building units around growing root')
xlabel('Distance from root (\mum)')
ylabel('Displacement normed by root expansion (%)')
imwrite(getframe(gcf).cdata, 'displacementEval.png')
%saveas(gcf, 'displacementEval.png')
%figure 
%geoBulk = reshape(b.bulkVector, [N,N]);
%imshow(geoBulk)

