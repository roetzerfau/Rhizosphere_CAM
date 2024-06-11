function [disp] = calculateDisplacementVectors(g,rootVector_ref, particleList, particleList_ref)
N = g.NX;
geoRoot = reshape(rootVector_ref, [N,N]);

D = bwdist(geoRoot);

%v = g.V0T(Ind,:);
%coord = mean(g.coordV(v,:),1);
%X = [coord; coordSource];
%d(i) = pdist(X,'euclidean');
disp = zeros(numel(particleList_ref),5);
for index = 1:numel(particleList_ref)
    particle_ref = particleList_ref{index};
    [coord_ref_row, coord_ref_col] = ind2sub([N,N],particle_ref);
    d = D(sub2ind([N,N], round(mean(coord_ref_row)), round(mean(coord_ref_col))));
    %d = mean(D(particle_ref));
    
    
%     v_ref = g.V0T(particle_ref,:);
%     coord_ref = g.coordV(v_ref,:);

    particle= particleList{index};
    [coord_row, coord_col] = ind2sub([N,N],particle);
%     v = g.V0T(particle,:);
%     coord = g.coordV(v,:);
    %% COM
    i_max = N;
    j_max = N;
    i = coord_ref_row;
    j = coord_ref_col;
    
    %% Row
    r_i = i_max/(2* pi);
    theta_i = (i/i_max) * 2 * pi; 
    x = r_i * cos(theta_i);
    y = j;
    z = r_i *sin(theta_i);
    
    
    
    x_bar = mean(x);
    y_bar = mean(y);
    z_bar = mean(z);
    
    theta_i = atan2(-z_bar,-x_bar) + pi;
    i_bar = i_max/ (2*pi) * theta_i;
    
    %% COl
    r_j = j_max/(2 * pi);
    theta_j = j/j_max * 2 * pi;
    x = i; 
    y = r_j * cos(theta_j);
    z = r_j * sin(theta_j);
    
    x_bar = mean(x);
    y_bar = mean(y);
    z_bar = mean(z);
    
    theta_j = atan2(-z_bar, - y_bar) + pi;
    j_bar = j_max/(2*pi) * theta_j;
    
    
    %% Displacement
    i_max = N;
    j_max = N;
    i_ref = coord_ref_row;
    i = coord_row;
    j_ref = coord_ref_col;
    j = coord_col;
    %% ROW
    r_i = i_max/(2* pi);
    theta_i = (i/i_max) * 2 * pi; 
    theta_i_ref = (i_ref/i_max) * 2 * pi; 
    
    x = r_i * cos(theta_i);
    y = j;
    z = r_i *sin(theta_i);
    
    x_ref = r_i * cos(theta_i_ref);
    y_ref = j_ref;
    z_ref = r_i *sin(theta_i_ref);
%     
%     x_d = x - x_ref;
%     y_d = y - y_ref;
%     z_d = z - z_ref;

    %theta_i = atan2(-z_d,-x_d) + pi;
    theta_i_ = atan2(-z,-x) + pi;
    theta_i_ref_ = atan2(-z_ref,-x_ref) + pi;
    
    
    theta_d1 = theta_i - theta_i_ref;
    theta_d2 = mod(theta_i - theta_i_ref - 2* pi, 2*pi);
    %theta_d2 = -sign(theta_d2).* mod(theta_d2, 2*pi);
    theta_d3 = theta_d2 - 2 * pi;
%     theta_d4 = theta_i_ref - theta_i - 2* pi;
 
    vectors = {theta_d1, theta_d2, theta_d3};
    theta_d = unique(cat(2, vectors{:}));
    %theta_d = union(theta_d1, theta_d2);
    i_d = i_max/ (2*pi)* theta_d;
    %% Col
    r_j = j_max/(2 * pi);
    theta_j = j/j_max * 2 * pi;
    theta_j_ref = j_ref/j_max * 2 * pi;
    
%     x = i; 
%     y = r_j * cos(theta_j);
%     z = r_j * sin(theta_j);
%     
%     x_ref = i_ref; 
%     y_ref = r_j * cos(theta_j_ref);
%     z_ref = r_j * sin(theta_j_ref);

    
%     theta_j = atan2(-z, - y) + pi;
%     theta_j_ref = atan2(-z_ref, - y_ref) + pi;

    theta_d1 = theta_j - theta_j_ref;
    theta_d2 = mod(theta_j - theta_j_ref - 2* pi, 2*pi);
    theta_d3 = theta_d2 - 2 * pi;
%     theta_d3 = theta_j_ref - theta_j;
%     theta_d4 = theta_j_ref - theta_j - 2 * pi;
    %theta_d = union(theta_d1, theta_d2);
    vectors = {theta_d1, theta_d2, theta_d3};
    theta_d = unique(cat(2, vectors{:}));
    j_d = j_max/(2*pi) * theta_d;
%% --
    %centOfGrav_row = mean(coord_row,1);
%     centOfGrav_ref_row =  mean(coord_ref_row);
%     centOfGrav_ref_col =  mean(coord_ref_col);
    %centOfGrav2 = coord_ref(1,:);
    
%     for j = 2:size(coord,1)
%         centOfGrav2 = (mod(centOfGrav2 +  coord_ref(j,:) +  N * N,N))/2;%/size(coord,1);
%     end
% 
%     centOfGrav2 = centOfGrav2
   
    
    
%     i
%     i_ref
%     i_d
    disp_x = coord_row - coord_ref_row;
    [B,I_i] = min(abs(i_d));
    
    disp_y = coord_col - coord_ref_col;
    [B,I_j] = min(abs(j_d));
%     un_i = unique(i_d);
%     un_j = unique(j_d);
    disp(index,1) = round(i_d(I_i));%disp_x(I);
    disp(index,2) = round(j_d(I_j));%disp_y(I);
    
    disp(index,3) = ceil(i_bar);
    disp(index,4) = ceil(j_bar);
     disp(index,5) = d;
    %disp(i,3) = d;
    %coord = mean(g.coordV(v,:),1);
end
