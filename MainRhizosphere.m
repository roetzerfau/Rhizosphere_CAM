function MainRhizosphere

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
numOuterIt  = 4000  ;    

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
rootCells_n_initial = 0;
rootCells_n_current = rootCells_n_initial;
rootCells_growingRate = 2;
isFirstRootCell = 1;

centerOfDomain = g.NX/2;
rootVector = 0 * ones(g.numT, 1); 
rootInd = centerOfDomain * g.NX + centerOfDomain;

rootInd_v = g.V0T(rootInd,:);
coord_root_ = g.coordV(rootInd_v,:);
coord_root = mean(coord_root_);

rootParticleList = {[rootInd]};
%rootVector(rootInd) = 1;
%in bulk vector alles drin
%bulkVector(rootInd) = 1;

stencil_root = stencil( g.NX, g.NX, rootInd, g.NX-2);%1 
vertices_ = (g.V0T(stencil_root,:));
vertices = reshape(transpose(vertices_),1,[]);
coord = g.coordV(vertices,:);
coord_t = coord;
coord = transpose(coord);
coord = reshape(coord,2,4,[]);
coord = mean(coord,2);
coord = transpose(coord(:,:));
%coord1 = reshape(coord,4,2,[]);
distance = sqrt(...
    (coord(:,1)-coord_root(1)).^2 +...
    (coord(:,2)-coord_root(2)).^2);



trapezNum = 1:g.numT;
vertices = reshape(transpose(g.V0T(trapezNum,:)),1,[]);
coord = transpose(g.coordV(vertices,:));
coord = reshape(coord,2,4,[]);
coord = mean(coord,2);
coord = transpose(coord(:,:));
distance = sqrt(...
    (coord(:,1)-coord_root(1)).^2 +...
    (coord(:,2)-coord_root(2)).^2);

[B,I] = sort(distance);
%erst mit Stencil nachbarn gucken und dann nochmal mit distance zu mitte
%sortieren
%% Creating Initial Distribution of Bulk and Biomass                        %% Can be changed
% visualizeDataEdges(g, edgeChargeVector, 'memoryEdges', 'edgeChargeVector', 0, 2);
% visualizeDataEdges(g, reactiveSurfaceVector, 'reactiveEdges', 'reactiveSurfaceVector', 0, 2);
% visualizeDataEdges(g, concPOMAgent, 'agent', 'concPOMAgent', 0, 2);
% visualizeDataEdges(g, POMagentAge, 'age', 'POMagentAge', 0, 2);
% visualizeDataSub(g, POMconcVector, 'POMconc', 'POMconc', 0);
% visualizeDataSub(g, POMageVector, 'POMage', 'POMage', 0); 
visualizeDataSub(g, bulkVector + POMVector + rootVector, 'cellType', 'solu', 0);
visualizeDataSub(g, rootVector, 'root', 'root', 0);
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

printInfoTUM(0,bulkVector,POMconcVector, concPOMAgent, edgeChargeVector, POMsolidEdgeList,...
    numFreePOMparticles, numEdgeTypes, totalPOMinputConc, totalPOMoutputConc, sumExcessPOM, POMagentInput,...
    POMocclusion_total, POMocclusion_attractive);

sumAgent = sum(concAgent);

for k = 1 : numOuterIt
%% Growing Root
rootCells_n_new = rootCells_growingRate * k;
%if(k < 0.66 * numOuterIt )
    rootCells_n_expected = rootCells_n_initial + rootCells_growingRate * k;
%else
  %  rootCells_n_expected = rootCells_n_expected - rootCells_growingRate;
%end
%if(rootCells_n_expected > rootCells_n_current)
    
% Mit Stencil        
     
%     for ind = 1:rootCells_n_expected
%          rootVector(stencil_root(ind)) = 1;
%           bulkVector(trapezNum(I(ind)))= 1;
%     end
% Mit Konzentrisch
% for ind = 1:rootCells_n_expected
%      rootVector(trapezNum(I(ind))) = 1;
%      bulkVector(trapezNum(I(ind)))= 1;
% end
% if(rootCells_n_expected < rootCells_n_current)
%      rootVector(trapezNum(I(rootCells_n_expected:rootCells_n_current))) = 0;
%      bulkVector(trapezNum(I(rootCells_n_expected:rootCells_n_current)))= 0;
% end
%rootCells_n_current = rootCells_n_expected;
%vielleicht noch mit Stencil gucken, ob Cellen direkt verbunden
if(k < 0.66 * numOuterIt )
    rootCell_index = 1;
    for ind = 1:rootCells_growingRate
        isnewRootCellFound = 0;
        while ~isnewRootCellFound 
            if(rootCell_index > g.NX^2)
                fprintf('No free Space found \n');
                break;
            end
            if(rootCell_index < 0)
                fprintf('Root disappeared \n');
                break;
            end
            entry = rootVector(trapezNum(I(rootCell_index)));
            entry_sten = stencil(g.NX,g.NX,trapezNum(I(rootCell_index)),1);
%             if(isFirstRootCell)
%                 isRootCellConnected = 1;
%             else
%                 s = sum(rootVector(entry_sten),'all');
%                 if(s > 0)
%                     isRootCellConnected = 1;
%                 else
%                     isRootCellConnected = 0;
%                 end
%             end
            
            if(entry == 0 && bulkVector(trapezNum(I(rootCell_index)))== 0 )%&& isRootCellConnected
                rootVector(trapezNum(I(rootCell_index))) = 1;
                bulkVector(trapezNum(I(rootCell_index))) = 1;
                isnewRootCellFound = 1;
            end
            rootCell_index = rootCell_index +1;
            if(isFirstRootCell && isnewRootCellFound == 1)
                isFirstRootCell = 0;
            end

        end
    end
else
    rootCell_index = size(I,1);
    for ind = 1:rootCells_growingRate
        islatestRootCellFound = 0;
        while ~islatestRootCellFound 
            if(rootVector(trapezNum(I(rootCell_index))) == 1)
                rootVector(trapezNum(I(rootCell_index))) = 0;
                bulkVector(trapezNum(I(rootCell_index))) = 0;
                islatestRootCellFound = 1;
            end
             rootCell_index = rootCell_index -1;
        end
       
    end
    
    
end
%vertices of Trap
%g.V0T
%coord of vertices 
%g.coordV
%end

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
        solidParticleList{ solidParticle },~] = moveParticles( particleSize, bigParticleStencilLayers_individual,...
        g, bulkVector, bulkTypeVector, particleTypeVector, POMVector, POMconcVector, POMageVector,...
        concAgent, concPOMAgent, POMagentAge, edgeChargeVector, reactiveSurfaceVector, ...
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
        POMParticleList{ POMParticle },~] = moveParticles( particleSize, bigParticleStencilLayers_individual,...
        g, bulkVector, bulkTypeVector,particleTypeVector, POMVector, POMconcVector, POMageVector,...
        concAgent, concPOMAgent, POMagentAge, edgeChargeVector, reactiveSurfaceVector,...
        NZd , fileID , POMParticleList{ POMParticle },sumAgent,4,0, attraction_type);  

end

% tic

%% Movement of aggregates
if false
% T_start = tic;
bigJumping = 1;
if bigJumping == 1
    
% identify all aggregates consisting of solid building units and POM
% particles
[particleList, particleContent] = particleInfoTUM(bulkVector, solidParticleList, POMParticleList);
for particle = 1 : length( particleList )
    particleSize = length( particleList{ particle } ); 
    if(size(particleContent{particle},1)<2)% kein Verbund
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
        concAgent, concPOMAgent, POMagentAge, edgeChargeVector, reactiveSurfaceVector, changedList,~] = ...
        moveParticles( particleSize, bigParticleStencilLayers_individual, g, bulkVector, bulkTypeVector, ... 
        particleTypeVector, POMVector, POMconcVector, POMageVector, concAgent, ...
        concPOMAgent, POMagentAge, edgeChargeVector, reactiveSurfaceVector, NZd , ...
        fileID ,particleList{ particle },sumAgent,4,0, attraction_type);  
  

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
    elseif plot_frequency == 1 && (k <= 5 || mod(k,100) == 0 || k == numOuterIt)
% elseif plot_frequency == 1 
%     uLagr       = projectDG2LagrangeSub( uDG );
%     visualizeDataSub(g, uLagr, 'u', 'solu', k);
    visualizeDataSub(g, bulkVector + POMVector + + rootVector, 'cellType', 'solu', k);
%     visualizeDataSub(g, POMconcVector, 'POMconc', 'POMconc', k);
%     visualizeDataSub(g, POMageVector, 'POMage', 'POMage', k);  
% %     visualizeDataEdges(g, concAgent, 'conc', 'agent', k);
%     visualizeDataEdges(g, edgeChargeVector, 'memoryEdges', 'edgeChargeVector', k, 2);
%     visualizeDataEdges(g, reactiveSurfaceVector, 'reactiveEdges', 'reactiveSurfaceVector', k, 2);
%     visualizeDataEdges(g, concPOMAgent, 'agent', 'concPOMAgent', k, 2);
%     visualizeDataEdges(g, POMagentAge, 'age', 'POMagentAge', k, 2);
% visualizeDataSub(g, particleTypeVector, 'particleType', 'solu', k); 
     %visualizeDataSub(g, bulkVector, 'bulkVector', 'solu', k);
     visualizeDataSub(g, rootVector, 'root', 'root', k);
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
