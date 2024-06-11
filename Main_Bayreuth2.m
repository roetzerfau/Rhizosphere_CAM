function Main_Bayreuth2
clc; clear;
%% Compiling the C++ Components, necessary to determine the particle size distribution
% mex particleSizeDistribution.cpp;

%% Definition of parameters               
plot_frequency = 1;  % 0: only initial and final state; 1: specified below

attraction_type = 5; % 1: old volume charges, 2: no charges, 3: edge Charges, 4: for TUM, 5: Freising paper

% Input files
inputMat = 'Input/BlankDomain_250.mat'; % contains initial state testMain250.mat   example_20.mat BlankDomain_20.mat  PaperConfigs/por05_34_500.mat  
inputMat = 'Input/C_gradient_50.mat';
%inputMat = 'Input/config.90.mat';
randomPOMinputShapes = 'Input/POMshapes250_15.mat'; % contains shapes of POM particles

% inputPOMmat = 'Input/POMinputTest.mat';
% can be used to give very specific POM input, e.g. which was given as an
% output from different simulation


inputTimeSteps = 'Input/inputParticleNum_125.mat';
% contains 'inputVector' with info how many POM particle should be added at
% every time step, randomly chosen from shapes given by
% 'randomPOMinputShapes'

inputRoot = 'Input/rootConfig.90.mat';
% Number of Time Steps
numOuterIt  = 1000;    

% Flag if POM decay should be considered (0: no, 1: yes)
POMdecayFlag = 1;

% Parameters
parameters.POMdecayRate = 0.0096; % 0.0351 day^-1: Recous1995 after 3 days;
% 0.0096 day^-1: Bucka 2019, 0.0239 day^-1: Recous1995 after 18 days;
parameters.startConcPOM = 1.4; % initial concentration in POM cells
parameters.POMminConc = 0.01; % minimal threshold for POM cells; below this
% value POM cells turn into pore cells
parameters.POMmaxConc = 1; % not relevant for now, but could be used for 
% spreading of excess POM (or bio phase) to neighboring cells above this
% value
parameters.carbonUseEfficiency = 0.45; % CENTURY model: 0.45 for below-ground,
% 0.55 for surface organic matter decomposition

parameters.POMagentDecayRate = 1.5 * parameters.POMdecayRate; % decomposition
% rate for agent
parameters.POMagentMin = 0.01; %CAUTION: concPOMagent can be below threshold,
% but memory edges are only created if conc. above threshold, and age when
% below threshold
parameters.POMagentMax = 1;

parameters.POMinputAfterNsteps = 10; % POM input is given after every
% N steps (high : 2, low : 10, no: 2000)
parameters.POMinputNumParticles = 1; % Number of POM particles given as 
% input after N steps (1 for 250 domain, 4 for 500 domain)
parameters.relocateFreePOMafterNsteps = 2000; % POM particles that were
% without attractive neighbor for N consecutive steps are relocated


parameters.rootGrowingRate =  0.0921;
parameters.rootShrinkingRate = 0.011;%ein drittel übrig


parameters.minConcMucilage = 0.1;
parameters.mucilageGrowingRate = 0.05;
parameters.mucilageDecayRate = 3.1;%7.9398;%0.22;%2.3026 ->sum/0.22
parameters.normalMucilageConcentration = 6.5;%-> immer kleiner als 0.005 *1.3 in rhizisphre
criticialValue = parameters.normalMucilageConcentration * 1.5;
parameters.constantMucilageDeposition = 333.8;
extraConcAmount = 0;
TrootGrowingBegin = 0;
TrootGrowingEnd = 100;%80
TrootShrinkingBegin = 100;
TrootShrinkingEnd = 1300;
TmucilageGrowingBegin = 0;
TmucilageGrowingEnd = 100;%80
k_start = 0;
procentRootMucilage = 1;
% add parameters concerning aging
% add parameters for stencil sizes
% add parameters for probabilities of breaking up


% POM particles are removed, if they have been free for more than
% 'removePOMthreshold' consecutive steps

removePOMthreshold = 1001;
sumAgent = 0;

%Flags:
isRoot = false
moveParticles = true


%% Creating domain for Simulation
[g, bulkVector,bulkTypeVector, concAgent, edgeChargeVector, reactiveSurfaceVector, particleTypeVector,...
    POMVector, POMconcVector, POMageVector, concPOMAgent, POMagentAge, solidParticleList, POMParticleList, randomPOMparticles, ...
     randomPOMparticlesSizes] = initializeDomainFromMat(inputMat, randomPOMinputShapes);

NZd = g.NX; 
totalPOMinputConc = sum(POMconcVector);
fileID = fopen( 'Move_bulk_log_file' , 'w' );
% fileID_1 = fopen( 'Verteilungen', 'w' );


%% Creating Initial Root
if false
    load(inputRoot,'rootVector', 'mucilageVector', 'mucilageConcVector', 'mucilageSurfaceVector' , 'MucilageagentAge','rootComplexList', 'rootComplexGraph', 'mucilageGraph');
else
    [rootVector, mucilageVector, mucilageConcVector,mucilageSurfaceVector, MucilageagentAge, rootComplexList, rootComplexGraph,mucilageGraph] = ...
creatingInitialRootComplex(g, bulkVector,isRoot);
bulkVector = bulkVector + rootVector;

end


%% Test
% concentrationVector = zeros(g.numT, 1);
% concentrationVector(1) =1;
% B = stencil( g.NX , g.NX , 1 , 4);
% A = stencil( g.NX , g.NX , 1 , 5);
% C = setdiff(A,B);
% 
% bulkVector(C) = 1;
% visualizeDataSub(g, concentrationVector, 'concentrationVector', 'concentrationVector', 0);
% concentrationVector = spreadConcentration(g, concentrationVector, bulkVector,1);
% visualizeDataSub(g, concentrationVector, 'concentrationVector', 'concentrationVector', 1);
%% Test End
%% Test 2
parameters.v_Cliquid = 1.2 * 10^-4;
parameters.K_Cliquid = 5 * 10^-4;
parameters.Resp_GE = 0.26;
parameters.Resp_Maint = 2.31*10^-6;%0.008;
parameters.eta_PAR = 1;


parameters.minConC_B = 0.0132;
parameters.initConC_B = 0.0539;
parameters.maxConcC_B = 0.3168;


parameters.C_N_S = 0.1;
parameters.C_N_B = 10;
parameters.C_N_POM = 100;
parameters.C_N_Root = 100;
parameters.N_initialMicrobes = 15;

parameters.tau_ode = 3600; %12 * 60;
C_SVector = ones(g.numT, 1) .* ~bulkVector * parameters.startConcPOM/100;%0.1
N_SVector =  C_SVector ./ parameters.C_N_S;

previousConcentration = sum(C_SVector)
Concentration_Phase = find(C_SVector > 0);
Concentration_Phase_Occupied = Concentration_Phase(bulkVector(Concentration_Phase) > 0 );
Concentration_Phase_NotOccupied = Concentration_Phase(bulkVector(Concentration_Phase) == 0);


C_BVector = zeros(g.numT, 1);
%centerOfDomain = g.NX/2;
%MBIntialCellInd = [centerOfDomain * g.NX + centerOfDomain, round(centerOfDomain/ 2) * g.NX + round(centerOfDomain),(centerOfDomain + round(centerOfDomain/ 2)) * g.NX + centerOfDomain ];
candidates = find(bulkVector == 0);
msize = numel(candidates);
MBIntialCellInd = candidates(randperm(msize, parameters.N_initialMicrobes));
%MBIntialCellInd = [];
C_BVector(MBIntialCellInd) = parameters.minConC_B;
N_BVector = C_BVector ./ parameters.C_N_B;

MB_Vector = C_BVector >= parameters.minConC_B;

MNVector = zeros(g.numT, 1);
% visualizeDataSub(g, C_BVector, 'C_BVector', 'C_BVector', 0);
% visualizeDataSub(g, C_SVector, 'C_SVector', 'C_SVector', 0);
% visualizeDataSub(g, N_BVector, 'N_BVector', 'N_BVector', 0);
% visualizeDataSub(g, N_SVector, 'N_SVector', 'N_SVector', 0);
parameters.mucilageGrowing = 1;
outerRootBorderInd = 1:250; 

for j = 1:100
    if(j > 25)
        parameters.mucilageGrowing = 0;
    end
T_C_N = tic; 
% [MB_Vector, N_SVector, C_SVector, N_BVector, C_BVector , MNVector, POMVector, POMconcVector] = ...
%     calculate_C_N(g, parameters, bulkVector, MB_Vector, N_SVector, C_SVector, N_BVector, C_BVector , MNVector, POMVector, POMconcVector, reactiveSurfaceVector, POMParticleList, POMageVector,outerRootBorderInd);
% fprintf('Time for C_N: %d \n', toc(T_C_N))
% numel(find(C_BVector > 0))
% C_BVector(find(C_BVector > 0))
% visualizeDataSub(g, C_BVector , 'C_BVector', 'C_BVector', j);
% visualizeDataSub(g, C_SVector * 1000, 'C_SVector', 'C_SVector', j);
% visualizeDataSub(g, N_BVector, 'N_BVector', 'N_BVector', j);
% visualizeDataSub(g, N_SVector * 100, 'N_SVector', 'N_SVector', j);

end
%% Test 2 End




global notConnectedEdgesValue;
notConnectedEdgesValue = g.NX* g.NX * 2;
pressureDistributionVector = zeros(g.numT, 1);
rootPressureDistributionVector = zeros(g.numT, 1);
mucilagePressureDistributionVector = zeros(g.numT, 1);

mantles = [0,10,20,30,40,50, 60 ,70];
porosity_table = zeros(numOuterIt, numel(mantles)-1);
%Obstacle
%bulkVector(N * N/4: N* (N/4+1)) = 1;
%bulkVector(N * N/2 + N/2 +3:N * N/2 + N/2 +10) = 1;

 
visualizeDataSub(g, bulkVector + POMVector + rootVector*2 + MB_Vector *4 , 'cellType', 'solu', k_start);
visualizeDataSub(g, C_BVector, 'C_BVector', 'C_BVector', 0);
visualizeDataSub(g, C_SVector, 'C_SVector', 'C_SVector', 0);
visualizeDataSub(g, N_BVector, 'N_BVector', 'N_BVector', 0);
visualizeDataSub(g, N_SVector, 'N_SVector', 'N_SVector', 0);

for k = k_start + 1 : numOuterIt
fprintf('k %d \n', k)

%% Root
currentAmountRootCells = sum(rootVector, 'all');
newAmountRootCells = 0;
if(TrootGrowingBegin < k && k <= TrootGrowingEnd && isRoot)
    newAmountRootCells = 900*(k- TrootGrowingBegin);%900*k;%100 * k;875
elseif(TrootShrinkingBegin < k && currentAmountRootCells > 1 &&isRoot)
    newAmountRootCells = floor(currentAmountRootCells * exp(-parameters.rootShrinkingRate));
else
    newAmountRootCells = currentAmountRootCells;
end
diffAmountRootCells = newAmountRootCells - currentAmountRootCells;
   
amountChangeCells = diffAmountRootCells;

rootPressureDistributionVector(:) = 0;
rootComplexList_old = rootComplexList;
rootEdges_old = g.CE0T(rootComplexList_old,:);
mucilageSurfaceVector(reshape(rootEdges_old,1,[])) = 0;
if(amountChangeCells >  0)
    T_growing = tic;
    [rootComplexGraph, bulkVector, rootComplexList, rootPressureDistributionVector...
     bulkTypeVector, particleTypeVector, POMVector, POMconcVector, POMageVector, ...
    concAgent, concPOMAgent, POMagentAge,MucilageagentAge, edgeChargeVector, reactiveSurfaceVector,mucilageSurfaceVector, mucilageVector,...
    POMParticleList, solidParticleList, newCellsInd] = ...
    rootMucilageComplexGrowing_3(g, rootComplexGraph, bulkVector, rootVector, rootComplexList, amountChangeCells,...
    bulkTypeVector, particleTypeVector, POMVector, POMconcVector, POMageVector, ...
    concAgent, concPOMAgent, POMagentAge,MucilageagentAge, edgeChargeVector, reactiveSurfaceVector, mucilageSurfaceVector,mucilageVector,...
    POMParticleList, solidParticleList,...
    attraction_type);
    fprintf('Time for RootGrowing: %d \n', toc(T_growing))
    
    rootVector(:) = 0;
    rootVector(rootComplexList) = 1;
elseif(amountChangeCells < 0)
    amountChangeCells = abs(amountChangeCells);
    deadCells = rootComplexList(end-amountChangeCells:end);
    [rootComplexGraph, bulkVector, rootVector, rootComplexList] = rootMucilageComplexShrinking(rootComplexGraph, bulkVector, rootVector, rootComplexList, deadCells);
    %vielleicht doch kein springen in mucialge
else

end
% if(parameters.constantMucilageDeposition > 0 && isRoot)
% 	notrootComplexList_l = cell(1,1);
% 	notrootComplexList_l{1} = find(rootVector == 0);
% 	mucilageRootEdgeList = calculatePOMsolidEdgeList(g, rootVector, ~rootVector, notrootComplexList_l);
% 	mucilageSurfaceVector(mucilageRootEdgeList{1}) = 1;
% 	if(k > TrootGrowingEnd + 50)
% 		randHelper = randi(500, size(mucilageRootEdgeList{1}));
% 		randHelper = randHelper <= 499;
% 		procentRootMucilage = procentRootMucilage * sum(randHelper)/numel(randHelper);
% 		nonMucilageEdge = floor(numel(randHelper) * (1-procentRootMucilage));
% 		mucilageSurfaceVector(mucilageRootEdgeList{1}(1:nonMucilageEdge)) = 0;
% 	end
% end



T_mucPress = tic;

highConcInd = sum(mucilageConcVector > criticialValue, 'all');
mucilagePressureDistributionVector(:) = 0;
if(numel(highConcInd)>0)
pressPoints  = find(bwdist(mucilageConcVector > criticialValue) == 1);
pressPoints = pressPoints(bulkVector(pressPoints) == 1);
mucilagePressureDistributionVector = calculatePressureDistribution(g, pressPoints, bulkVector, rootComplexList);
end
M = cat(2,rootPressureDistributionVector(:),mucilagePressureDistributionVector(:));
pressureDistributionVector = reshape(max(M,[],2),size(rootPressureDistributionVector));
fprintf('Time for mucPressure: %d \n', toc(T_mucPress))


%-----------------------------------
%wegschieben, aber auch höhere Attraktivität an möglichen zukünftigen
%Randpunkten/mucilagee andockpunkte kennzeichnen


%% Doing the POM decay
if POMdecayFlag == 1

T_spread = tic;

%C_SVector = easyDiffusiveStep(g, C_SVector, bulkVector, 5);

%[bulkVector, bulkTypeVector, POMVector, POMconcVector, POMageVector, POMParticleList, C_SVector, N_SVector] = calculateonlyPOMdecay(g, parameters, bulkVector, bulkTypeVector, POMVector, POMconcVector,...
 %    reactiveSurfaceVector, POMParticleList, POMageVector, C_SVector, N_SVector);



fprintf('Time for POM decay and spreading of agent: %d \n', toc(T_spread))

%% Input of POM particles
% input of POM particles using specified parameters

% if mod(k,parameters.POMinputAfterNsteps) == 0
%     occupiedCells = bulkVector + mucilageVector;
%     for inputParticle = 1 : parameters.POMinputNumParticles
%     [occupiedCells, bulkTypeVector, POMVector, POMconcVector, POMageVector, POMParticleList, totalPOMinputConc] = placePOMparticleRandomly(g, occupiedCells, bulkTypeVector, ...
%         POMVector, POMconcVector, POMageVector, POMParticleList, randomPOMparticles, randomPOMparticlesSizes, totalPOMinputConc);
%     end
%     bulkVector = occupiedCells - mucilageVector;
% end


% aging of memory edges
% POMagentAge = calculatePOMagentAge(parameters, POMagentAge, edgeChargeVector, concPOMAgent );
% [edgeChargeVector, POMagentAge] = randomAging(POMagentAge, edgeChargeVector);
% 
% MucilageagentAge = MucilageagentAge + 1;
% [mucilageSurfaceVector, MucilageagentAge] = randomAging(MucilageagentAge, mucilageSurfaceVector);

end


%if moveParticles

%% Movement of solid building units
T_start = tic;
% move inseparable solid building units
for solidParticle = 1 : length( solidParticleList )
    particleSize = length( solidParticleList{ solidParticle } ); 

    % calculate stencil of solid building units depending on its area
    % (=particleSize)
%     bigParticleStencilLayers_individual = 1;
    bigParticleStencilLayers_individual = min(5, ceil(20/(particleSize)^0.5));
    meanPressure = mean(pressureDistributionVector(solidParticleList{ solidParticle }));
    bigParticleStencilLayers_individual = max(ceil(meanPressure * 5), bigParticleStencilLayers_individual);
    %bulkMucilageVector = bulkVector + mucilageVector;
%   apply CAM for single solid building unit
    [bulkVector,bulkTypeVector, particleTypeVector, POMVector, POMconcVector, POMageVector,...
        concAgent, concPOMAgent, POMagentAge, MucilageagentAge, edgeChargeVector, reactiveSurfaceVector,mucilageSurfaceVector,...
        mucilageVector, pressureDistributionVector,...
        solidParticleList{ solidParticle },~] = moveParticles( particleSize, bigParticleStencilLayers_individual,...
        g, bulkVector, bulkTypeVector, particleTypeVector, POMVector, POMconcVector, POMageVector,...
        concAgent, concPOMAgent, POMagentAge, MucilageagentAge, edgeChargeVector, reactiveSurfaceVector,mucilageSurfaceVector, ...
        mucilageVector, pressureDistributionVector, ...
        NZd , fileID ,solidParticleList{ solidParticle },sumAgent,4,0, attraction_type,0);  
     %bulkVector = bulkMucilageVector - mucilageVector;
end

%% Movement of POM particles
% move all POM particles
for POMParticle = 1 : length( POMParticleList )
    particleSize = length( POMParticleList{ POMParticle } );  
    
    % calculate stencil of POM particle depending on its area
%     bigParticleStencilLayers_individual = 1;
    bigParticleStencilLayers_individual = min(5, ceil(20/(particleSize)^0.5));
     meanPressure = mean(pressureDistributionVector(POMParticleList{ POMParticle }));
    bigParticleStencilLayers_individual = max(ceil(meanPressure * 5), bigParticleStencilLayers_individual);
    % apply CAM for single POM particle
    bulkMBVector = bulkVector + MB_Vector;
    [bulkMBVector,bulkTypeVector, particleTypeVector, POMVector, POMconcVector, POMageVector,...
        concAgent, concPOMAgent, POMagentAge, MucilageagentAge, edgeChargeVector, reactiveSurfaceVector, mucilageSurfaceVector,...
        mucilageVector, pressureDistributionVector,...
        POMParticleList{ POMParticle },~] = moveParticles( particleSize, bigParticleStencilLayers_individual,...
        g, bulkMBVector, bulkTypeVector,particleTypeVector, POMVector, POMconcVector, POMageVector,...
        concAgent, concPOMAgent, POMagentAge,MucilageagentAge, edgeChargeVector, reactiveSurfaceVector,mucilageSurfaceVector,...
        mucilageVector, pressureDistributionVector, ...
        NZd , fileID , POMParticleList{ POMParticle },sumAgent,4,0, attraction_type,0);  
     bulkVector = bulkMBVector - MB_Vector;
end

% tic

%% Movement of aggregates
if true
% T_start = tic;
bigJumping = 1;
if bigJumping == 1
    
% identify all aggregates consisting of solid building units and POM
% particles
[particleList, particleContent] = particleInfoTUM(bulkVector-(rootVector), solidParticleList, POMParticleList);
for particle = 1 : length( particleList )
    particleSize = length( particleList{ particle } ); 
    if(size(particleContent{particle},1)<2)% kein Verbund
        continue
    end
    if(size(particleContent{particle},1) == 10)% contains root
        continue
    end
    test_ind = particleList{particle}(1);
    
    bigParticleStencilLayers_individual = min(5, ceil(20/(particleSize)^0.5));
    meanPressure = mean(rootPressureDistributionVector( particleList{ particle }));
    bigParticleStencilLayers_individual = max(ceil(meanPressure * 5), bigParticleStencilLayers_individual);
    
    % bigParticleStencilLayers_individual = 1;
    if particleSize > 20000
        bigParticleStencilLayers_individual = 0;
    end
    if bigParticleStencilLayers_individual == 0
        continue
    end
    
     bulkMBVector = bulkVector + MB_Vector;
     
    [bulkMBVector,bulkTypeVector, particleTypeVector, POMVector, POMconcVector, POMageVector, ...
        concAgent, concPOMAgent, POMagentAge,MucilageagentAge, edgeChargeVector, reactiveSurfaceVector,mucilageSurfaceVector,...
        mucilageVector, pressureDistributionVector,...
        changedList,~] = ...
        moveParticles( particleSize, bigParticleStencilLayers_individual, g, bulkMBVector, bulkTypeVector, ... 
        particleTypeVector, POMVector, POMconcVector, POMageVector, concAgent, ...
        concPOMAgent, POMagentAge, MucilageagentAge, edgeChargeVector, reactiveSurfaceVector, mucilageSurfaceVector, ...
        mucilageVector, pressureDistributionVector, ...    
        NZd , fileID ,particleList{ particle },sumAgent,4,0, attraction_type,0);  
    bulkVector = bulkMBVector - MB_Vector;

    if test_ind ~= changedList(1) % Listen müssen angepasst werden
        sten_ind = find(stencil(g.NX,g.NX,test_ind,bigParticleStencilLayers_individual)==changedList(1));
        for part = 1: size(particleContent{particle},1)  
            switch particleContent{particle}(part,1)
                case 1 
                    sten = stencil(g.NX,g.NX,solidParticleList{ particleContent{particle}(part,2)},bigParticleStencilLayers_individual);
                    solidParticleList{particleContent{particle}(part,2)}=sten(:,sten_ind)';
                case 2 
                    sten = stencil(g.NX,g.NX,POMParticleList{ particleContent{particle}(part,2)},bigParticleStencilLayers_individual);
                    POMParticleList{particleContent{particle}(part,2)}=sten(:,sten_ind)';
            end
        end
    end 
  
end
end


fprintf('Time for relocation: %d \n', toc(T_start))
end
%end



%% C_N
if(k > 25)
        parameters.mucilageGrowing = 0;
else
        parameters.mucilageGrowing = 1;
end

T_C_N = tic; 
[MB_Vector, N_SVector, C_SVector, N_BVector, C_BVector , MNVector, POMVector, POMconcVector] = ...
    calculate_C_N(g, parameters, bulkVector, MB_Vector, N_SVector, C_SVector, N_BVector, C_BVector , MNVector, POMVector, POMconcVector, reactiveSurfaceVector, POMParticleList, POMageVector,outerRootBorderInd);
fprintf('Time for C_N: %d \n', toc(T_C_N))
numel(find(C_BVector > 0))
C_BVector(find(C_BVector > 0))

%% Porosity
% 
% [porosity_t] = calculatePorosity(g,mantles,rootVector, bulkVector - rootVector);
% porosity_table(k,:) = porosity_t';
% writematrix(porosity_table,'FinalConfig/porosity_table.xls')
%% Displacement Vectors


%%
T2 = tic;


    if plot_frequency == 1 && (k <= 200 || mod(k,25) == 0 || k == numOuterIt)
    visualizeDataSub(g, bulkVector + POMVector + rootVector*2 + MB_Vector *4 , 'cellType', 'solu', k);
    visualizeDataSub(g, C_BVector , 'C_BVector', 'C_BVector', k);
    visualizeDataSub(g, C_SVector, 'C_SVector', 'C_SVector', k);
    visualizeDataSub(g, N_BVector, 'N_BVector', 'N_BVector', k);
    visualizeDataSub(g, N_SVector, 'N_SVector', 'N_SVector', k);
    end


%%
%     if( mod(k,5) == 0 || k == numOuterIt)
% %         particleListHelper = particleList;
%         particleList = solidParticleList;
%         fileName    = ['FinalConfig/config','.', num2str(k),'.mat']; 
%         save(fileName,'g','bulkVector','bulkTypeVector','POMconcVector', 'concPOMAgent','edgeChargeVector','POMagentAge',...
%             'POMVector', 'POMageVector', 'POMParticleList', 'particleList', 'reactiveSurfaceVector', 'particleTypeVector')        
% %         particleList = particleListHelper;
% 
%         fileName    = ['FinalConfig/rootConfig','.', num2str(k),'.mat']; 
%         save(fileName,'rootVector', 'mucilageVector', 'mucilageConcVector','mucilageSurfaceVector', 'MucilageagentAge', 'rootComplexList', 'rootComplexGraph','mucilageGraph')       
%     end

    
end  % for k
fclose( fileID );
% fclose( fileID_1);
end
