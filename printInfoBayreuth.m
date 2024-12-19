function printInfoBayreuth(k, N_SVector, C_SVector, N_BVector, C_BVector ,C_MNVector, N_MNVector, C_POMconcVector,N_POMconcVector,CO2Vector, leakedNVector,...
    CUE, f_C,f_BD_C,B_C,R, f_N,f_BD_N,B_N, sumleakedN_S, f_POM_C, f_POM_N, f_MN_C, f_MN_N ,POMageVector,soil_particle_NNZ, folder_output)

if k == 0
   flag = 'w';
else
   flag = 'a';
end

fileName    =  folder_output + '/C_BVector.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, sum(C_BVector));
fclose(fileID_k);

fileName    = folder_output + '/C_BVectorNNZ.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, sum(C_BVector > 0));
fclose(fileID_k);



fileName    =  folder_output + '/C_SVector.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, sum(C_SVector));
fclose(fileID_k);

fileName    =  folder_output + '/C_SVectorNNZ.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, sum(C_SVector > 0));
fclose(fileID_k);




fileName    =  folder_output + '/C_POMconcVector.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, sum(C_POMconcVector));
fclose(fileID_k);

fileName    =  folder_output + '/C_POMconcVectorNNZ.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, sum(C_POMconcVector > 0));
fclose(fileID_k);



fileName    =  folder_output + '/C_MNVector.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, sum(C_MNVector));
fclose(fileID_k);

fileName    =  folder_output + '/C_MNVectorNNZ.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, sum(C_MNVector > 0));
fclose(fileID_k);



fileName    =  folder_output + '/CO2Vector.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, sum(CO2Vector));
fclose(fileID_k);

fileName    =  folder_output + '/CO2VectorNNZ.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, sum(CO2Vector > 0));
fclose(fileID_k);



total_C = sum(C_BVector) + sum(C_SVector) + sum(C_POMconcVector) + sum(C_MNVector)+sum(CO2Vector);
fileName    =  folder_output + '/total_C.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, total_C);
fclose(fileID_k);

total_C = sum((C_BVector + C_SVector + C_POMconcVector + C_MNVector +CO2Vector) > 0);
fileName    =  folder_output + '/total_CNNZ.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, total_C);
fclose(fileID_k);



fileName    = folder_output + '/N_BVector.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, sum(N_BVector));
fclose(fileID_k);

fileName    = folder_output + '/N_BVectorNNZ.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, sum(N_BVector > 0));
fclose(fileID_k);


fileName    =  folder_output + '/N_SVector.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, sum(N_SVector));
fclose(fileID_k);

fileName    =  folder_output + '/N_SVectorNNZ.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, sum(N_SVector > 0));
fclose(fileID_k);

fileName    =  folder_output + '/N_POMconcVector.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, sum(N_POMconcVector));
fclose(fileID_k);

fileName    =  folder_output + '/N_POMconcVectorNNZ.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, sum(N_POMconcVector > 0));
fclose(fileID_k);


fileName    = folder_output + '/N_MNVector.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, sum(N_MNVector));
fclose(fileID_k);

fileName    = folder_output + '/N_MNVectorNNZ.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, sum(N_MNVector > 0));
fclose(fileID_k);



fileName    =  folder_output + '/leakedNVector.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, sum(leakedNVector));
fclose(fileID_k);

fileName    =  folder_output + '/leakedNVectorNNZ.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, sum(leakedNVector > 0));
fclose(fileID_k);


total_N = sum(N_BVector) + sum(N_SVector) + sum(N_POMconcVector) + sum(N_MNVector) + sumleakedN_S;
fileName    =  folder_output + '/total_N.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, total_N);
fclose(fileID_k);

total_N = sum(N_BVector) + sum(N_SVector) + sum(N_POMconcVector) + sum(N_MNVector) + sumleakedN_S;
fileName    =  folder_output + '/total_NNNZ.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, total_N);
fclose(fileID_k);


fileName    =  folder_output + '/CUE.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, CUE);
fclose(fileID_k);

fileName    =  folder_output + '/R.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, R);
fclose(fileID_k);

fileName    =  folder_output + '/f_B_C.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, B_C);
fclose(fileID_k);

fileName    =  folder_output + '/f_BD_C.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, f_BD_C);
fclose(fileID_k);

fileName    =  folder_output + '/f_C.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, f_C);
fclose(fileID_k);



fileName = folder_output + '/sumleakedN_S.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, sumleakedN_S);
fclose(fileID_k);

fileName    =  folder_output + '/f_B_N.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, B_N);
fclose(fileID_k);

fileName    =  folder_output + '/f_BD_N.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, f_BD_N);
fclose(fileID_k);

fileName    =  folder_output + '/f_N.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, f_N);
fclose(fileID_k);



fileName = folder_output + '/f_POM_C.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, f_POM_C);
fclose(fileID_k);

fileName    =  folder_output + '/f_POM_N.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, f_POM_N);
fclose(fileID_k);

fileName    =  folder_output + '/f_MN_C.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, f_MN_C);
fclose(fileID_k);

fileName    =  folder_output + '/f_MN_N.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, f_MN_N);
fclose(fileID_k);





fileName    =  folder_output + '/soil_particleNNZ.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, soil_particle_NNZ);
fclose(fileID_k);

indices_new = find(C_MNVector > 0 & POMageVector ~= max(POMageVector));
indices_old = find(C_MNVector > 0 & POMageVector == max(POMageVector));
printequal = sum(C_MNVector) - (sum(C_MNVector(indices_new)) + sum(C_MNVector(indices_old)))

fileName    =  folder_output + '/C_MNVector_new.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, sum(C_MNVector(indices_new)));
fclose(fileID_k);

fileName    =  folder_output + '/C_MNVector_newNNZ.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, sum(C_MNVector(indices_new) > 0));
fclose(fileID_k);

fileName    =  folder_output + '/C_MNVector_old.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, sum(C_MNVector(indices_old)));
fclose(fileID_k);

fileName    =  folder_output + '/C_MNVector_oldNNZ.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, sum(C_MNVector(indices_old) > 0));
fclose(fileID_k);


indices_new = find(N_MNVector > 0 & POMageVector ~= max(POMageVector));
indices_old = find(N_MNVector > 0 & POMageVector == max(POMageVector));
if(sum(C_MNVector) - (sum(C_MNVector(indices_new)) + sum(C_MNVector(indices_old)))> 0.000001)
    fprintf("falsch print %f \n ", sum(C_MNVector) - (sum(C_MNVector(indices_new)) + sum(C_MNVector(indices_old))))
end

fileName    =  folder_output + '/N_MNVector_new.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, sum(N_MNVector(indices_new)));
fclose(fileID_k);

fileName    =  folder_output + '/N_MNVector_newNNZ.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, sum(N_MNVector(indices_new) > 0));
fclose(fileID_k);

fileName    =  folder_output + '/N_MNVector_old.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, sum(N_MNVector(indices_old)));
fclose(fileID_k);

fileName    =  folder_output + '/N_MNVector_oldNNZ.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, sum(N_MNVector(indices_old) > 0));
fclose(fileID_k);

end

