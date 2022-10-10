function EigenValueOnPorosity(  )

%% Compiling the C++ Components, necessary to determine the particle size distribution
mex particleSizeDistribution.cpp;

%% Defining Porosity Vector and Parameters
POR         = 0.5 * ones(1,10);
NX          = 25;
numSmallPor = 0.05;
numBigPor   = 0.05;

%% Fixed Parameters
NZu         = 0;
width       = NX;
intBound    = NX * ones(NX+1,1);
upBound     = NX * ones(NX+1,1);

%% Computation
g = createDomainFolded(width, NX, NX, NZu, intBound, upBound);

lambda_1    = zeros(size(POR));
lambda_2    = zeros(size(POR));

i = 1;

while i <= length(POR)
    porosity = POR(i);
    helperPor = POR(i);
    while true
        bulkVector = createBulkVector(g, floor(NX*NX*helperPor*numBigPor), floor(NX*NX*helperPor*numSmallPor));
        if sum(~bulkVector) < NX^2 * (porosity - 0.01)
            helperPor = helperPor + 0.01;
        elseif sum(~bulkVector) > NX^2 * (porosity + 0.01)
           helperPor = helperPor - 0.01;
        else
            break
        end
    end    
    K = SimpleHomogenizationFunc( bulkVector );
    if sum(sum(isnan(K))) > 0
        continue;
    end
    a   = eigs(K);
    if a(1) > 1.1 || a(2)>1.1
        continue;
    end
    lambda_1(i)   = a(1);
    lambda_2(i)   = a(2);
    
    particleSizeD = particleSizeDistribution(bulkVector.');
    
    fileName    = ['particleSizes.', num2str(i)];
    file        = fopen(fileName, 'wt');
    fprintf(file,'%d %d\n', particleSizeD.');
    fclose(file);
    
    i = i + 1;

end

fid = fopen('lambda', 'wt' );
fprintf( fid, '%d %18.14f %18.14f\n', [(1:length(POR)).', lambda_1.', lambda_2.'].' );
fclose( fid );

end