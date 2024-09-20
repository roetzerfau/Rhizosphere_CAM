function printInfoBayreuth(k, N_SVector, C_SVector, N_BVector, C_BVector ,C_MNVector, N_MNVector, C_POMconcVector,N_POMconcVector,CO2Vector, sumleakedN_S,CUE, folder_output)

if k == 0
   flag = 'w';
else
   flag = 'a';
end

fileName    =  folder_output + '/C_BVector.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, sum(C_BVector));
fclose(fileID_k);

fileName    =  folder_output + '/C_SVector.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, sum(C_SVector));
fclose(fileID_k);

fileName    =  folder_output + '/C_POMconcVector.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, sum(C_POMconcVector));
fclose(fileID_k);

fileName    =  folder_output + '/C_MNVector.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, sum(C_MNVector));
fclose(fileID_k);

fileName    =  folder_output + '/CO2Vector.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, sum(CO2Vector));
fclose(fileID_k);

total_C = sum(C_BVector) + sum(C_SVector) + sum(C_POMconcVector) + sum(C_MNVector)+sum(CO2Vector);
fileName    =  folder_output + '/total_C.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, total_C);
fclose(fileID_k);

fileName    = folder_output + '/N_BVector.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, sum(N_BVector));
fclose(fileID_k);

fileName    =  folder_output + '/N_SVector.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, sum(N_SVector));
fclose(fileID_k);

fileName    =  folder_output + '/N_POMconcVector.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, sum(N_POMconcVector));
fclose(fileID_k);

fileName    = folder_output + '/N_MNVector.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, sum(N_MNVector));
fclose(fileID_k);

fileName    =  folder_output + '/sumleakedN_S.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, sumleakedN_S);
fclose(fileID_k);

total_N = sum(N_BVector) + sum(N_SVector) + sum(N_POMconcVector) + sum(N_MNVector) + sumleakedN_S;
fileName    =  folder_output + '/total_N.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, total_N);
fclose(fileID_k);


fileName    =  folder_output + '/CUE.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, CUE);
fclose(fileID_k);
end

