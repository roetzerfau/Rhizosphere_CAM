function printInfoBayreuth(k, N_SVector, C_SVector, N_BVector, C_BVector ,C_MNVector, N_MNVector, C_POMconcVector,N_POMconcVector,  name)

if k == 0
   flag = 'w';
else
   flag = 'a';
end

fileName    =  'txtdata'+ string(name)+ '/C_BVector.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, sum(C_BVector));
fclose(fileID_k);

fileName    =  'txtdata'+ string(name)+ '/C_SVector.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, sum(C_SVector));
fclose(fileID_k);

fileName    =  'txtdata'+ string(name)+ '/C_POMconcVector.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, sum(C_POMconcVector));
fclose(fileID_k);

fileName    =  'txtdata'+ string(name)+ '/C_MNVector.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, sum(C_MNVector));
fclose(fileID_k);

total_C = sum(C_BVector) + sum(C_SVector) + sum(C_POMconcVector) + sum(C_MNVector);
fileName    =  'txtdata'+ string(name)+ '/total_C.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, total_C);
fclose(fileID_k);

fileName    =  'txtdata'+ string(name)+ '/N_BVector.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, sum(N_BVector));
fclose(fileID_k);

fileName    =  'txtdata'+ string(name)+ '/N_SVector.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, sum(N_SVector));
fclose(fileID_k);

fileName    =  'txtdata'+ string(name)+ '/N_POMconcVector.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, sum(N_POMconcVector));
fclose(fileID_k);

fileName    =  'txtdata'+ string(name)+ '/N_MNVector.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, sum(N_MNVector));
fclose(fileID_k);

total_N = sum(N_BVector) + sum(N_SVector) + sum(N_POMconcVector) + sum(N_MNVector);
fileName    =  'txtdata'+ string(name)+ '/total_N.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, total_N);
fclose(fileID_k);
end

