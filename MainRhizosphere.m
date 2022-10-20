function MainRhizosphere
clc; clear;
%% Compiling the C++ Components, necessary to determine the particle size distribution
% mex particleSizeDistribution.cpp;

%% Definition of parameters               
plot_frequency = 1;  % 0: only initial and final state; 1: specified below

attraction_type = 5; % 1: old volume charges, 2: no charges, 3: edge Charges, 4: for TUM, 5: Freising paper

% Input files
inputMat = 'Input/testMain250.mat'; % contains initial state
randomPOMinputShapes = 'Input/POMshapes250_15.mat'; % contains shapes of POM particles

% inputPOMmat = 'Input/POMinputTest.mat';
% can be used to give very specific POM input, e.g. which was given as an
% output from different simulation

inputTimeSteps = 'Input/inputParticleNum_125.mat';
% contains 'inputVector' with info how many POM particle should be added at
% every time step, randomly chosen from shapes given by
% 'randomPOMinputShapes'

% Number of Time Steps
numOuterIt  = 1000  ;    

% Flag if POM decay should be considered (0: no, 1: yes)
POMdecayFlag = 1;

% Parameters
parameters.POMdecayRate = 0.0096; % 0.0351 day^-1: Recous1995 after 3 days;
% 0.0096 day^-1: Bucka 2019, 0.0239 day^-1: Recous1995 after 18 days;
parameters.startConcPOM = 1; % initial concentration in POM cells
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

parameters.POMinputAfterNsteps = 2000; % POM input is given after every
% N steps (high : 2, low : 10, no: 2000)
parameters.POMinputNumParticles = 1; % Number of POM particles given as 
% input after N steps (1 for 250 domain, 4 for 500 domain)
parameters.relocateFreePOMafterNsteps = 2000; % POM particles that were
% without attractive neighbor for N consecutive steps are relocated

% add parameters concerning aging
% add parameters for stencil sizes
% add parameters for probabilities of breaking up


% POM particles are removed, if they have been free for more than
% 'removePOMthreshold' consecutive steps
removePOMthreshold = 1001;


%% Creating domain for Simulation
[g, bulkVector,bulkTypeVector, concAgent, edgeChargeVector, reactiveSurfaceVector, particleTypeVector,...
    POMVector, POMconcVector, POMageVector, concPOMAgent, POMagentAge, solidParticleList, POMParticleList, randomPOMparticles, ...
     randomPOMparticlesSizes] = initializeDomainFromMat(inputMat, randomPOMinputShapes);
%bulkVector enthält POM und solidParticles
% load(inputPOMmat, 'inputPOMtime', 'inputPOMparticles', 'inputPOMparticlesConc'); %loads: inputPOMtime, inputPOMparticles, inputPOMparticlesConc

load(inputTimeSteps, 'inputVector'); %loads inputVector
 
% geo = zeros(10,10);
% porosity = 0.5; 
%  [g, bulkVector,bulkTypeVector, concAgent, edgeChargeVector, reactiveSurfaceVector, particleTypeVector,...
%      POMVector, POMconcVector, POMageVector, concPOMAgent, POMagentAge, solidParticleList, POMParticleList]...
%      = initializeDomainSimple(geo, porosity);  

% list that contains for every POM particle the indices of neighboring
% solid edes
POMsolidEdgeList = calculatePOMsolidEdgeList(g, bulkVector, POMVector, POMParticleList);

totalPOMinputConc = sum(POMconcVector);
totalPOMoutputConc = 0;
sumExcessPOM = 0;
POMagentInput = 0;

removedPOMparticles = {};
timeRemovedPOMparticles = [];
removedPOMparticlesConc = {};

NZd = g.NX; 

fileID = fopen( 'Move_bulk_log_file' , 'w' );
% fileID_1 = fopen( 'Verteilungen', 'w' );


%% Creating Initial Root
rootNr = 1;
rootCells_n_initial = 0;
rootCellsCurrentAmount = rootCells_n_initial;
rootCells_growingRate = 10;
rootCells_shrinkingRate = -10;
rootCellsExpectedAmount = rootCellsCurrentAmount + rootCells_growingRate;
isRootGrowing = true;

rootVector = 0 * ones(g.numT, 1);
rootParticleList = {[]};
rootPressureEdgeVector = zeros(g.numCE,1);

N = g.NX;
notConnectedEdgesValue = N*N;

% [I,J]=ndgrid(1:N,1:N);
% IJ=[I(:),J(:)];
% D=pdist2(IJ,IJ);
% A=(D>=0.001 & D<=sqrt(2)*1.001).*D .*notConnectedEdgesValue; %adjacency matrix

diagVec1 = sparse(repmat([ones(N-1, 1); 0], N, 1));  % Make the first diagonal vector
                                             %   (for horizontal connections)
diagVec1 = diagVec1(1:end-1);                % Remove the last value
diagVec2 = sparse([0; diagVec1(1:(N*(N-1)))] * sqrt(2));       % Make the second diagonal vector
                                             %   (for anti-diagonal connections)
diagVec3 = sparse(ones(N*(N-1), 1));                 % Make the third diagonal vector
                                             %   (for vertical connections)
diagVec4 = sparse(diagVec2(2:end-1));                % Make the fourth diagonal vector
                                             %   (for diagonal connections)
adj = diag(diagVec1, 1)+...                  % Add the diagonals to a zero matrix
      diag(diagVec2, N-1)+...
      diag(diagVec3, N)+...
      diag(diagVec4, N+1);
adj = adj+adj.';                             % Add the matrix to a transposed copy of
                                             %   itself to make it symmetric

adj = adj .*notConnectedEdgesValue;

%ff = sum(adj == A,'all')

rootGraph = graph(adj);
clear adj diagVec1 diagVec2 diagVec3 diagVec4
 
%find nearest pore space to center of domain 
centerOfDomain = g.NX/2;
rootIntialCellInd = centerOfDomain * g.NX + centerOfDomain;
%d = distances(rootGraph);
%d = d(rootIntialCellInd,:);
[TR,d] = shortestpathtree(rootGraph,rootIntialCellInd);
[sortedd, I] = sort(d);
j =find(bulkVector(I) == 0,1);
rootNewCellsInd = I(j);
rootIntialCellInd = I(j);
rootDeadCellsInd = [];


indices = 1:(N*N);
indices = reshape(indices,[N N]);

%Obstacle
%bulkVector(N * N/4: N* (N/4+1)) = 1;
%bulkVector(N * N/2 + N/2 +3:N * N/2 + N/2 +10) = 1;
%% Creating Initial Distribution of Bulk and Biomass                        %% Can be changed
% visualizeDataEdges(g, edgeChargeVector, 'memoryEdges', 'edgeChargeVector', 0, 2);
% visualizeDataEdges(g, reactiveSurfaceVector, 'reactiveEdges', 'reactiveSurfaceVector', 0, 2);
% visualizeDataEdges(g, concPOMAgent, 'agent', 'concPOMAgent', 0, 2);
% visualizeDataEdges(g, POMagentAge, 'age', 'POMagentAge', 0, 2);
visualizeDataEdges(g, rootPressureEdgeVector, 'pressureEdges', 'rootPressureEdgeVector', 0,2);
% visualizeDataSub(g, POMconcVector, 'POMconc', 'POMconc', 0);
% visualizeDataSub(g, POMageVector, 'POMage', 'POMage', 0); 
visualizeDataSub(g, bulkVector + POMVector + rootVector, 'cellType', 'solu', 0);
%visualizeDataSub(g, rootVector, 'root', 'root', 0);
numEdgeTypes =  countEdgeTypes(g, bulkVector, POMVector, solidParticleList, ...
    edgeChargeVector, reactiveSurfaceVector, particleTypeVector);
numEdgeTypesPOMparticles =  countEdgeTypesPOMparticles(g, bulkVector, POMVector,...
    POMParticleList, edgeChargeVector, reactiveSurfaceVector);

POMocclusion_total = 1 - sum(numEdgeTypesPOMparticles(:,1))/sum(sum(numEdgeTypesPOMparticles));
POMocclusion_attractive = (sum(numEdgeTypesPOMparticles(:,4)) + ...
    sum(numEdgeTypesPOMparticles(:,5)))/sum(sum(numEdgeTypesPOMparticles));
occlusionFactor = (numEdgeTypesPOMparticles(:,4)+numEdgeTypesPOMparticles(:,5))./sum(numEdgeTypesPOMparticles,2);

%%% count number of free POM particles: Ignore for now
indFreePOMparticles = [];
    
for i = 1 : length(POMParticleList)
    if sum(edgeChargeVector(POMsolidEdgeList{i})) + sum(reactiveSurfaceVector(POMsolidEdgeList{i})) == 0
       indFreePOMparticles = [indFreePOMparticles; i];
    end
end

indFreePOMparticles = find(occlusionFactor < 0.05);

freePOMparticles = [indFreePOMparticles ones(length(indFreePOMparticles),1)];
freePOMparticlesOld = freePOMparticles;

numFreePOMparticles = length(indFreePOMparticles);
%%%

%printInfoTUM(0,bulkVector,POMconcVector, concPOMAgent, edgeChargeVector, POMsolidEdgeList,...
    %numFreePOMparticles, numEdgeTypes, totalPOMinputConc, totalPOMoutputConc, sumExcessPOM, POMagentInput,...
    %POMocclusion_total, POMocclusion_attractive);

sumAgent = sum(concAgent);

for k = 1 : numOuterIt
%% Growing Root

if isRootGrowing
    for i = 1:size(rootNewCellsInd)
        if(bulkVector(rootNewCellsInd(i)) == 0)
            neighInd = neighbors(rootGraph, rootNewCellsInd(i));    
            rootBorderCellsInd = intersect(neighInd,rootParticleList{rootNr});
            
            h = rootNewCellsInd(i) * ones(size(rootBorderCellsInd,1),1);
            edgeInd = findedge(rootGraph, rootBorderCellsInd, h);

            %nur van neumman weiterwachsen -> man dürfte nur Kanten mit 1 teilen
            rootGraph.Edges.Weight(edgeInd) = rootGraph.Edges.Weight(edgeInd)./notConnectedEdgesValue;
        
            rootVector(rootNewCellsInd(i)) = 1;
            bulkVector(rootNewCellsInd(i)) = 1;
            rootParticleList{rootNr} = ...
            [rootParticleList{rootNr} rootNewCellsInd(i)];
        else
            fprintf('Cell cant grow \n');
        end
        
    end    
else
    for i = 1:size(rootDeadCellsInd)
        
            neighInd = neighbors(rootGraph, rootDeadCellsInd(i));    
            rootNeighborCellsInd = neighInd(~ismember(neighInd,rootParticleList{rootNr}));
            
            h = rootDeadCellsInd(i) * ones(size(rootNeighborCellsInd,1),1);
            edgeInd = findedge(rootGraph, rootNeighborCellsInd, h);

            %nur van neumman weiterwachsen -> man dürfte nur Kanten mit 1 teilen
            rootGraph.Edges.Weight(edgeInd) = rootGraph.Edges.Weight(edgeInd).*notConnectedEdgesValue;
            
            rootVector(rootDeadCellsInd(i)) = 0;
            bulkVector(rootDeadCellsInd(i)) = 0;
            rootParticleList{rootNr} = ...
            rootParticleList{rootNr}(rootParticleList{rootNr} ~=rootDeadCellsInd(i));
        
    end
    
    
end
rootCellsCurrentAmount = size(rootParticleList{rootNr},2);
if(k < 0.66 * numOuterIt )
    rate = rootCells_growingRate;
else
    rate = rootCells_shrinkingRate;
end

rootCellsExpectedAmount = rootCellsExpectedAmount + rate;
rootGrowingPotential = rootCellsExpectedAmount - rootCellsCurrentAmount;
if(rootGrowingPotential>0)
    isRootGrowing = true;
    %-------------------------------
    rootFreeCellsInd = find(rootVector ~= 1);
    %sort the potential new rootCells by their possibility to grow
    %(distance to initial source cell)
    %d = distances(rootGraph);
    %d = d(rootFreeCellsInd,rootIntialCellInd);
    [TR,d] = shortestpathtree(rootGraph,rootFreeCellsInd,rootIntialCellInd);
    [sortedd, I] = sort(d);
    rootFreeCellsInd = rootFreeCellsInd(I);
    if(size(rootFreeCellsInd)<size(rootGrowingPotential,1))
        rootGrowingPotential = size(rootFreeCellsInd,1);
    end
    rootNewCellsInd = rootFreeCellsInd(1:rootGrowingPotential);

    %--------------------------------------
    st = stencil( g.NX, g.NX,rootNewCellsInd , 1);

    [rootSurfaceEdgeList] = getSolidSurfaceEdges(g,rootParticleList, rootVector);
    rootPressureEdgeVector(:) = 0;
    %rootPressureEdgeVector(rootSurfaceEdgeList{rootNr}) = 1;
    growingCellEdges = g.CE0T(st,:);
    pressEdgesInd= intersect(rootSurfaceEdgeList{rootNr},growingCellEdges );
    rootPressureEdgeVector(pressEdgesInd) = 1;

elseif(rootGrowingPotential<0)
    isRootGrowing = false;
    %-------------------------------
    %sort the potential new rootCells by their possibility to grow
    %(distance to initial source cell)
    rootcurrentCellsInd = rootParticleList{rootNr};
    %d = distances(rootGraph);
    %d = d(rootcurrentCellsInd,rootIntialCellInd);
    [TR,d] = shortestpathtree(rootGraph,rootcurrentCellsInd,rootIntialCellInd);
    [sortedd, I] = sort(d);
    rootcurrentCellsInd = rootcurrentCellsInd(I);
    rootDeadCellsInd = rootcurrentCellsInd(end:end +1 + rootGrowingPotential);
    rootPressureEdgeVector(:) = 0;
else
end
%----------------------------------------------
%alle cellen außer root
%distances mit rootGraph berechen
%alle RootCellen rauschmeißen
%sortieren
%dann 

%abfolge:
%erst versuchen sachen zu setzen, danach neue Nachbarn berechnen und
%PressureEdges machen
%beim ersten mal rootNeighbourCellsInd definitv auf etwas richtiges setzen
if true
%% Doing the POM decay
if POMdecayFlag == 1

T_spread = tic;
[bulkVector, bulkTypeVector, POMVector, POMconcVector, POMageVector, concPOMAgent, edgeChargeVector, POMParticleList, POMsolidEdgeList, sumExcessPOM, POMagentInput] = ...
    calculatePOMdecay(g, parameters, bulkVector, bulkTypeVector, POMVector, POMconcVector, POMageVector, concPOMAgent, edgeChargeVector,...
    reactiveSurfaceVector, POMParticleList, sumExcessPOM, POMagentInput);
fprintf('Time for POM decay and spreading of agent: %d \n', toc(T_spread))

%% Input of POM particles
% input of POM particles using specified parameters
if mod(k,parameters.POMinputAfterNsteps) == 0
    for inputParticle = 1 : parameters.POMinputNumParticles
    [bulkVector, bulkTypeVector, POMVector, POMconcVector, POMageVector, POMParticleList, totalPOMinputConc] = placePOMparticleRandomly(g, bulkVector, bulkTypeVector, ...
        POMVector, POMconcVector, POMageVector, POMParticleList, randomPOMparticles, randomPOMparticlesSizes, totalPOMinputConc);
    end
end

% input of POM particles using 'inputTimeSteps'
% for inputParticle = 1 : inputVector(k)
%     [bulkVector, bulkTypeVector, POMVector, POMconcVector, POMageVector, POMParticleList, totalPOMinputConc] = placePOMparticleRandomly(g, bulkVector, bulkTypeVector, ...
%         POMVector, POMconcVector, POMageVector, POMParticleList, randomPOMparticles, randomPOMparticlesSizes, totalPOMinputConc);
% end


% % relocate free POM particles
% if mod(k,parameters.relocateFreePOMafterNsteps) == 0
%  [bulkVector, bulkTypeVector, POMVector, POMconcVector, POMageVector, POMParticleList] = ...
%      replaceFreePOMparticles(g, bulkVector, bulkTypeVector, edgeChargeVector, ...
%      reactiveSurfaceVector, POMVector, POMconcVector, POMageVector, POMParticleList);
% end

% aging of memory edges
POMagentAge = calculatePOMagentAge(parameters, POMagentAge, edgeChargeVector, concPOMAgent );
[edgeChargeVector, POMagentAge] = randomAging(POMagentAge, edgeChargeVector);

% % remove free POM particles above threshold
% [bulkVector, bulkTypeVector, POMVector, POMconcVector, POMageVector, POMParticleList, ...
%      removedPOMparticles, removedPOMparticlesConc, timeRemovedPOMparticles, totalPOMoutputConc] = ...
%      removeFreePOMparticles(g, bulkVector, bulkTypeVector, edgeChargeVector, reactiveSurfaceVector, ...
%      POMVector, POMconcVector, POMageVector, POMParticleList, removePOMthreshold, freePOMparticles, ...
%      removedPOMparticles, removedPOMparticlesConc, timeRemovedPOMparticles, totalPOMoutputConc, k);
 
 
% % add input POM particles using 'inputPOMmat'
% indsInput = find(inputPOMtime == k);
% 
% for j = 1 : length(indsInput)%bulkVector
%     [bulkVector, bulkTypeVector, POMVector, POMconcVector, POMageVector, POMParticleList, totalPOMinputConc] = ...
%      placeNewPOMparticle(g, bulkVector, bulkTypeVector, POMVector, POMconcVector, POMageVector, POMParticleList, ...
%      inputPOMparticles{indsInput(j)}, 10, inputPOMparticlesConc{indsInput(j)}, totalPOMinputConc);
% end

end
   

%% Movement of solid building units
T_start = tic;
% move inseparable solid building units
for solidParticle = 1 : length( solidParticleList )
    particleSize = length( solidParticleList{ solidParticle } ); 

    % calculate stencil of solid building units depending on its area
    % (=particleSize)
%     bigParticleStencilLayers_individual = 1;
    bigParticleStencilLayers_individual = min(5, ceil(20/(particleSize)^0.5));
    
%   apply CAM for single solid building unit
    [bulkVector,bulkTypeVector, particleTypeVector, POMVector, POMconcVector, POMageVector,...
        concAgent, concPOMAgent, POMagentAge, edgeChargeVector, reactiveSurfaceVector,...
        rootVector, rootPressureEdgeVector,...
        solidParticleList{ solidParticle },~] = moveParticles( particleSize, bigParticleStencilLayers_individual,...
        g, bulkVector, bulkTypeVector, particleTypeVector, POMVector, POMconcVector, POMageVector,...
        concAgent, concPOMAgent, POMagentAge, edgeChargeVector, reactiveSurfaceVector, ...
        rootVector, rootPressureEdgeVector, ...
        NZd , fileID ,solidParticleList{ solidParticle },sumAgent,4,0, attraction_type);  

end

%% Movement of POM particles
% move all POM particles
for POMParticle = 1 : length( POMParticleList )
    particleSize = length( POMParticleList{ POMParticle } );  
    
    % calculate stencil of POM particle depending on its area
%     bigParticleStencilLayers_individual = 1;
    bigParticleStencilLayers_individual = min(5, ceil(20/(particleSize)^0.5));
    
    % apply CAM for single POM particle
    [bulkVector,bulkTypeVector, particleTypeVector, POMVector, POMconcVector, POMageVector,...
        concAgent, concPOMAgent, POMagentAge, edgeChargeVector, reactiveSurfaceVector, ...
        rootVector, rootPressureEdgeVector,...
        POMParticleList{ POMParticle },~] = moveParticles( particleSize, bigParticleStencilLayers_individual,...
        g, bulkVector, bulkTypeVector,particleTypeVector, POMVector, POMconcVector, POMageVector,...
        concAgent, concPOMAgent, POMagentAge, edgeChargeVector, reactiveSurfaceVector,...
        rootVector, rootPressureEdgeVector, ...
        NZd , fileID , POMParticleList{ POMParticle },sumAgent,4,0, attraction_type);  

end

% tic

%% Movement of aggregates
if true
% T_start = tic;
bigJumping = 1;
if bigJumping == 1
    
% identify all aggregates consisting of solid building units and POM
% particles
[particleList, particleContent] = particleInfoTUM(bulkVector-rootVector, solidParticleList, POMParticleList, rootParticleList);
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
    % bigParticleStencilLayers_individual = 1;
    if particleSize > 20000
        bigParticleStencilLayers_individual = 0;
    end
    if bigParticleStencilLayers_individual == 0
        continue
    end
    [bulkVector,bulkTypeVector, particleTypeVector, POMVector, POMconcVector, POMageVector, ...
        concAgent, concPOMAgent, POMagentAge, edgeChargeVector, reactiveSurfaceVector,...
        rootVector, rootPressureEdgeVector,...
        changedList,~] = ...
        moveParticles( particleSize, bigParticleStencilLayers_individual, g, bulkVector, bulkTypeVector, ... 
        particleTypeVector, POMVector, POMconcVector, POMageVector, concAgent, ...
        concPOMAgent, POMagentAge, edgeChargeVector, reactiveSurfaceVector,  ...
        rootVector, rootPressureEdgeVector, ...    
        NZd , fileID ,particleList{ particle },sumAgent,4,0, attraction_type);  
  

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
end
%%
T2 = tic;

% uDG(find(bulkVector == 0),1) = 0;
% uDG(find(bulkVector == 1),1) = 1;
% uDG(find((bulkVector == 1) & (particleTypeVector == 0)),1) = 2;

    %alt
    if plot_frequency == 0 && k == numOuterIt
%     uLagr       = projectDG2LagrangeSub( uDG );
    visualizeDataSub(g, bulkVector + POMVector + rootVector, 'cellType', 'solu', k);
    visualizeDataSub(g, POMconcVector, 'POMconc', 'POMconc', k);
    visualizeDataSub(g, POMageVector, 'POMage', 'POMage', k);  
%     visualizeDataEdges(g, concAgent, 'conc', 'agent', k);
    visualizeDataEdges(g, edgeChargeVector, 'memoryEdges', 'edgeChargeVector', k, 2);
    visualizeDataEdges(g, reactiveSurfaceVector, 'reactiveEdges', 'reactiveSurfaceVector', k, 2);
    visualizeDataEdges(g, concPOMAgent, 'agent', 'concPOMAgent', k, 2);
    visualizeDataEdges(g, POMagentAge, 'age', 'POMagentAge', k, 2);
% visualizeDataSub(g, particleTypeVector, 'particleType', 'solu', k); 
% visualizeDataSub(g, bulkVector, 'bulkVector', 'solu', k); 
    elseif plot_frequency == 1 && (k <= 20 || mod(k,10) == 0 || k == numOuterIt)
% elseif plot_frequency == 1 
%     uLagr       = projectDG2LagrangeSub( uDG );
%     visualizeDataSub(g, uLagr, 'u', 'solu', k);
    visualizeDataSub(g, bulkVector + POMVector + + rootVector, 'cellType', 'solu', k);
%     visualizeDataSub(g, POMconcVector, 'POMconc', 'POMconc', k);
%     visualizeDataSub(g, POMageVector, 'POMage', 'POMage', k);  
% %     visualizeDataEdges(g, concAgent, 'conc', 'agent', k);
%     visualizeDataEdges(g, edgeChargeVector, 'memoryEdges', 'edgeChargeVector', k, 2);
%     visualizeDataEdges(g, reactiveSurfaceVector, 'reactiveEdges', 'reactiveSurfaceVector', k, 2);
    visualizeDataEdges(g, rootPressureEdgeVector, 'pressureEdges', 'rootPressureEdgeVector', k,2);
%     visualizeDataEdges(g, concPOMAgent, 'agent', 'concPOMAgent', k, 2);
%     visualizeDataEdges(g, POMagentAge, 'age', 'POMagentAge', k, 2);
% visualizeDataSub(g, particleTypeVector, 'particleType', 'solu', k); 
     %visualizeDataSub(g, bulkVector, 'bulkVector', 'solu', k);
     %visualizeDataSub(g, rootVector, 'root', 'root', k);
    end

%% Some posprocessing: Ignore for now
numEdgeTypes =  countEdgeTypes(g, bulkVector, POMVector, solidParticleList, ...
    edgeChargeVector, reactiveSurfaceVector, particleTypeVector);
numEdgeTypesPOMparticles =  countEdgeTypesPOMparticles(g, bulkVector, POMVector,...
    POMParticleList, edgeChargeVector, reactiveSurfaceVector);

POMocclusion_total = 1 - sum(numEdgeTypesPOMparticles(:,1))/sum(sum(numEdgeTypesPOMparticles));
POMocclusion_attractive = (sum(numEdgeTypesPOMparticles(:,4)) + ...
    sum(numEdgeTypesPOMparticles(:,5)))/sum(sum(numEdgeTypesPOMparticles));
occlusionFactor = (numEdgeTypesPOMparticles(:,4)+numEdgeTypesPOMparticles(:,5))./sum(numEdgeTypesPOMparticles,2);

POMsolidEdgeList = calculatePOMsolidEdgeList(g, bulkVector, POMVector, POMParticleList);

%%%
freePOMparticlesOld = freePOMparticles;

indFreePOMparticles = [];
    
% for i = 1 : length(POMParticleList)
%     if sum(edgeChargeVector(POMsolidEdgeList{i})) + sum(reactiveSurfaceVector(POMsolidEdgeList{i})) == 0
%        indFreePOMparticles = [indFreePOMparticles; i];
%     end
% end

indFreePOMparticles = find(occlusionFactor < 0.07);

freePOMparticles = [indFreePOMparticles zeros(length(indFreePOMparticles),1)];

[Lia, Locb] = ismember(freePOMparticles(:,1), freePOMparticlesOld(:,1));
freePOMparticles(find(Lia),:) = freePOMparticlesOld(Locb(find(Locb)), :);
if(~isempty(freePOMparticles))
    freePOMparticles(:,2) = freePOMparticles(:,2) + 1;
end

numFreePOMparticles = length(indFreePOMparticles);
%%

% printInfoTUM(k,bulkVector,POMconcVector, concPOMAgent, edgeChargeVector, POMsolidEdgeList,...
%     numFreePOMparticles, numEdgeTypes, totalPOMinputConc, totalPOMoutputConc, sumExcessPOM, POMagentInput,...
%     POMocclusion_total, POMocclusion_attractive);
    if(k < 50 || mod(k,50) == 0 || k == numOuterIt)
%         particleListHelper = particleList;
        particleList = solidParticleList;
        fileName    = ['FinalConfig/config','.', num2str(k),'.mat']; 
        save(fileName,'g','bulkVector','bulkTypeVector','POMconcVector', 'concPOMAgent','edgeChargeVector','POMagentAge',...
            'POMVector', 'POMageVector', 'POMParticleList', 'particleList', 'reactiveSurfaceVector', 'particleTypeVector',...
            'removedPOMparticles', 'removedPOMparticlesConc', 'timeRemovedPOMparticles')        
%         particleList = particleListHelper;
    end

    
end  % for k
fclose( fileID );
% fclose( fileID_1);
end
