function printInfoBayreuth(k,C_BVector,C_SVector, MNVector, name)

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

fileName    =  'txtdata'+ string(name)+ '/MNVector.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, sum(MNVector));
fclose(fileID_k);

end

