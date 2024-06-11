fileName = 'adsf';
tLvl = 3;
fileName    = ['vtk/' , fileName, '.', num2str(tLvl), '.vtu'];
fileID = fopen(fileName ,'wt');
fprintf(fileID,"asdf")
fclose(fileID)