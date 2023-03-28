clear; close all;
a = load('config.5.mat');
b = load('config.100.mat');

b_root = load('rootConfig.100.mat');
a_root = load('rootConfig.5.mat');
c =  b.particleList(1:2);
disp = calculateDisplacementVectors(a.g, a_root.rootVector, b.particleList(1:100), a.particleList(1:100));
N = a.g.NX;
geoBulk = reshape(a.bulkVector, [N,N]);
%geoBulk = cat(3, geoBulk, geoBulk, geoBulk);
figure


 
X = disp(:,3);
Y = disp(:,4);
U = disp(:,1);
V = disp(:,2);
%V(:) = 0;
d = disp(:,5);
ind = sub2ind(size(geoBulk), X, Y);
%ind = a.particleList{1};
[row,col] = ind2sub(size(geoBulk),ind);
% row = X(1);
% col = Y(2);
%geoBulk(row,col) = 0.5;
% geoBulk(row,col,1) = 1;
% geoBulk(row,col,2:3) = 0;
geoBulk(ind) = 0.5;
imshow(geoBulk)
hold on
quiver(Y,X,V,U, 0)
figure
M = sqrt(U.^2 + V.^2);
scatter(d, M);
figure 
geoBulk = reshape(b.bulkVector, [N,N]);
imshow(geoBulk)
