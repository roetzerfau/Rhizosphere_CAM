function InitialStateTUM

%% Compiling the C++ Components, necessary to determine the particle size distribution
% mex particleSizeDistribution.cpp;

%% Definition of parameters                %% Can be changed
bilder = 1;                                % 0: nur Bild von Anfangs- und Endzustand, 1: Bild zu jedem Schritt, 2: Bild zu ersten 10 Zuständen und Endzustand, 3: Bild nach allen 100 Schritten

attraction_type = 5; % 1: old volume charges, 2: no charges, 3: edge Charges, 4: for TUM, 5: Freising paper

randomPOMinputTxt = 'Input/POMshapes250_15.mat';
inputMat = 'Input/testMain250.mat';

% Number of Time Steps
numOuterIt  = 500;     

% Parameters
parameters.POMdecayRate = 0.;
parameters.startConcPOM = 1;
parameters.POMminConc = 0;
parameters.POMmaxConc = 1;
parameters.carbonUseEfficiency = 0.1;
parameters.POMagentDecayRate = 0.01;
parameters.POMagentMin = 0;
parameters.POMagentMax = 1;
parameters.POMinputAfterNsteps = 2000;
parameters.POMinputNumParticles = 0;
parameters.relocateFreePOMafterNsteps = 2000;


%% Creating domain for Simulation
[g, bulkVector,bulkTypeVector, concAgent, edgeChargeVector, reactiveSurfaceVector, particleTypeVector,...
    POMVector, POMconcVector, POMageVector, concPOMAgent, POMagentAge, solidParticleList, POMParticleList, randomPOMparticles, ...
     randomPOMparticlesSizes] = initializeDomainFromMat(inputMat, randomPOMinputTxt);
 POMsolidEdgeList = calculatePOMsolidEdgeList(g, bulkVector, POMVector, POMParticleList);
% [POMParticleList, ~] = createPOMParticleList(POMVector);
% solidParticleList =  createSolidParticleList (particleTypeVector);
totalPOMinputConc = sum(POMconcVector);
sumExcessPOM = 0;
POMagentInput = 0;

NZd = g.NX; 
markEbdr        = g.idE == 8 | g.idE == 8 | g.idE == 8 | g.idE == 8 | g.idE == 8;

fileID = fopen( 'Move_bulk_log_file' , 'w' );
%% Creating Initial Distribution of Bulk and Biomass                        %% Can be changed
visualizeDataEdges(g, edgeChargeVector, 'memoryEdges', 'edgeChargeVector', 0, 2);
visualizeDataEdges(g, reactiveSurfaceVector, 'reactiveEdges', 'reactiveSurfaceVector', 0, 2);
visualizeDataEdges(g, concPOMAgent, 'agent', 'concPOMAgent', 0, 2);
visualizeDataEdges(g, POMagentAge, 'age', 'POMagentAge', 0, 2);
visualizeDataSub(g, POMconcVector, 'POMconc', 'POMconc', 0);  
visualizeDataSub(g, POMageVector, 'POMage', 'POMage', 0);
visualizeDataSub(g, bulkVector + POMVector, 'particleType', 'particleType', 0); 

numEdgeTypes =  countEdgeTypes(g, bulkVector, POMVector, solidParticleList, ...
    edgeChargeVector, reactiveSurfaceVector, particleTypeVector);
printInfoTUM(0,bulkVector,POMconcVector, concPOMAgent, edgeChargeVector, POMsolidEdgeList,...
    numEdgeTypes, totalPOMinputConc, sumExcessPOM, POMagentInput);
printInfoTUM(k,bulkVector,POMconcVector, concPOMAgent, edgeChargeVector, POMsolidEdgeList,...
    numFreePOMparticles, numEdgeTypes, totalPOMinputConc, totalPOMoutputConc, sumExcessPOM, POMagentInput,...
    POMocclusion_total, POMocclusion_attractive)

%% Calculation of outer Loop for Movement of Bulk and Biomass Growth/Shrinkage

for k = 1 : numOuterIt
    fprintf('Step %d of %d\n', k, numOuterIt)
%% Add new POM particles
if mod(k,parameters.POMinputAfterNsteps) == 0
    for inputParticle = 1 : parameters.POMinputNumParticles
    [bulkVector, bulkTypeVector, POMVector, POMconcVector, POMageVector, POMParticleList, totalPOMinputConc] = placePOMparticleRandomly(g, bulkVector, bulkTypeVector, ...
        POMVector, POMconcVector, POMageVector, POMParticleList, randomPOMparticles, randomPOMparticlesSizes, totalPOMinputConc);
%     [bulkVector, bulkTypeVector, POMVector, POMconcVector, POMParticleList] = placePOMparticleRandomly(g, bulkVector, bulkTypeVector, ...
%         POMVector, POMconcVector, POMParticleList, randomPOMparticles, randomPOMparticlesSizes); 
    end
%     [POMParticleList, ~] = createPOMParticleList(POMVector);
end
if mod(k,parameters.relocateFreePOMafterNsteps) == 0
 [bulkVector, bulkTypeVector, POMVector, POMconcVector, POMageVector, POMParticleList] = ...
     replaceFreePOMparticles(g, bulkVector, bulkTypeVector, edgeChargeVector, ...
     reactiveSurfaceVector, POMVector, POMconcVector, POMageVector, POMParticleList);
end    
     
% Update POM age
    if k > 1
        POMageVector(POMageVector > 0) = POMageVector(POMageVector > 0) + 1;
        POMageVector(POMVector == 0) = 0;
    end
    

%% Doing the movement of Bulk
assert(sumAgent == sum(concAgent),'sumAgent')

% assert(numSolid == sum(bulkVector), 'solid Bloecke stimmen nicht')
T_start = tic;
% move solid particles

for solidParticle = 1 : length( solidParticleList )
    particleSize = length( solidParticleList{ solidParticle } ); 

%     bigParticleStencilLayers_individual = 1;
    bigParticleStencilLayers_individual = min(5, ceil(20/(particleSize)^0.5));
%     solidParticle
    assert(sumAgent==sum(concAgent),'sumAgent')
    [bulkVector,bulkTypeVector, particleTypeVector, POMVector, POMconcVector, POMageVector, uDG, q1DG, q2DG, uDGbac, q1DGbac,...
        q2DGbac, concAgent, concPOMAgent, POMagentAge, edgeChargeVector, reactiveSurfaceVector,...
        solidParticleList{ solidParticle },~] = moveParticles( particleSize, bigParticleStencilLayers_individual,...
        g, bulkVector, bulkTypeVector, particleTypeVector, POMVector, POMconcVector, POMageVector, uDG, q1DG, q2DG, uDGbac,...
        q1DGbac, q2DGbac, concAgent, concPOMAgent, POMagentAge, edgeChargeVector, reactiveSurfaceVector, 0, ...
        markEbdr, N, bioMvector, NZd , fileID ,solidParticleList{ solidParticle },sumAgent,4,0, attraction_type);  
    assert(sumAgent==sum(concAgent),'sumAgent')

end

%% move POM particles


for POMParticle = 1 : length( POMParticleList )
    particleSize = length( POMParticleList{ POMParticle } ); 

%     bigParticleStencilLayers_individual = 1;
    bigParticleStencilLayers_individual = min(5, ceil(20/(particleSize)^0.5));
    
    assert(sumAgent==sum(concAgent),'sumAgent')
    [bulkVector,bulkTypeVector, particleTypeVector, POMVector, POMconcVector, POMageVector, uDG, q1DG, q2DG, uDGbac, q1DGbac,...
        q2DGbac, concAgent, concPOMAgent, POMagentAge, edgeChargeVector, reactiveSurfaceVector, ...
        POMParticleList{ POMParticle },~] = moveParticles( particleSize, 5, g, bulkVector, bulkTypeVector,...
        particleTypeVector, POMVector, POMconcVector, POMageVector, uDG, q1DG, q2DG, uDGbac, q1DGbac, q2DGbac, concAgent, ...
        concPOMAgent, POMagentAge, edgeChargeVector, reactiveSurfaceVector, 0, markEbdr, N, bioMvector, NZd , ...
        fileID , POMParticleList{ POMParticle },sumAgent,4,0, attraction_type);  
    assert(sumAgent==sum(concAgent),'sumAgent')

end

% tic

% move big particles, typeflag = 4
T_start = tic;
bigJumping = 1;
if bigJumping == 1
[particleList, particleContent] = particleInfoTUM(bulkVector, solidParticleList, POMParticleList);
for particle = 1 : length( particleList )
    particleSize = length( particleList{ particle } ); 
    if(size(particleContent{particle},1)<2)% kein Verbund
        continue
    end
    test_ind = particleList{particle}(1);
    
%     bigParticleStencilLayers_individual = round(bigParticleStencilLayersConstant*1/(particleSize)^0.5);
    bigParticleStencilLayers_individual = min(5, ceil(20/(particleSize)^0.5));
% bigParticleStencilLayers_individual = 1;
    if particleSize > 20000
        bigParticleStencilLayers_individual = 0;
    end
    if bigParticleStencilLayers_individual == 0
        continue
    end
    assert(sumAgent==sum(concAgent),'sumAgent')
    [bulkVector,bulkTypeVector, particleTypeVector, POMVector, POMconcVector, POMageVector, uDG, q1DG, q2DG, uDGbac, q1DGbac, ...
        q2DGbac, concAgent, concPOMAgent, POMagentAge, edgeChargeVector, reactiveSurfaceVector, changedList,~] = ...
        moveParticles( particleSize, bigParticleStencilLayers_individual, g, bulkVector, bulkTypeVector, ... 
        particleTypeVector, POMVector, POMconcVector, POMageVector, uDG, q1DG, q2DG, uDGbac, q1DGbac, q2DGbac, concAgent, ...
        concPOMAgent, POMagentAge, edgeChargeVector, reactiveSurfaceVector, 0, markEbdr, N, bioMvector, NZd , ...
        fileID ,particleList{ particle },sumAgent,4,0, attraction_type);  
    assert(sumAgent==sum(concAgent),'sumAgent')
  
    % switch inds in QuartzLits, GoethitList, IlliteList
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

    assert(sumAgent==sum(concAgent),'sumAgent_rotation_before')
end
end

fprintf('Time for jumping of aggregates: %d ', toc(T_start))
size(find(particleTypeVector))
size(find(bulkVector))
%%
T2 = tic;
% solidParticleList = createSolidParticleList (particleTypeVector);
% fprintf('Time for createSolidParticleLists: %d ', toc(T2))

%     uDG(bulkTypeVector == QuartzCharge,1) = -1;
%     uDG(bulkVector == 1,1) = -1;
    
%     uDG(bulkTypeVector == -1,1) = -1;
%     uDG(bulkTypeVector == 1,1) = 0;
    uDG(find(bulkVector == 0),1) = 0;
    uDG(find(bulkVector == 1),1) = 1;
    uDG(find((bulkVector == 1) & (particleTypeVector == 0)),1) = 2;
%     size(find(bulkVector == 1))
%     size(find((bulkVector == 1) & (particleTypeVector == 0)))
    %alt
    if bilder == 0 && k == numOuterIt
%     uLagr       = projectDG2LagrangeSub( uDG );
%     visualizeDataSub(g, uLagr, 'u', 'solu', k);
    visualizeDataSub(g, particleTypeVector, 'particleType', 'solu', k);
    elseif bilder == 1 && (k <= 5 || mod(k,100) == 0 || k == numOuterIt)
% elseif bilder == 1 
    uLagr       = projectDG2LagrangeSub( uDG );
    visualizeDataSub(g, uLagr, 'u', 'solu', k);
    visualizeDataSub(g, POMconcVector, 'POMconc', 'POMconc', k);  
    visualizeDataSub(g, POMageVector, 'POMage', 'POMage', k);  
%     visualizeDataEdges(g, concAgent, 'conc', 'agent', k);
visualizeDataEdges(g, edgeChargeVector, 'conc', 'edgeChargeVector', k, 2);
visualizeDataEdges(g, reactiveSurfaceVector, 'conc', 'reactiveSurfaceVector', k, 2);
visualizeDataEdges(g, concPOMAgent, 'agent', 'concPOMAgent', k, 2);
visualizeDataEdges(g, POMagentAge, 'age', 'POMagentAge', k, 2);
% visualizeDataSub(g, particleTypeVector, 'particleType', 'solu', k); 
% visualizeDataSub(g, bulkVector, 'bulkVector', 'solu', k); 
    elseif bilder == 2 && (k <= 10 || k == numOuterIt)	
%     uLagr       = projectDG2LagrangeSub( uDG );
%     visualizeDataSub(g, uLagr, 'u', 'solu', k); 
    visualizeDataSub(g, particleTypeVector, 'particleType', 'solu', k);  
    elseif bilder == 3 && (mod(k,100) == 0 || k == numOuterIt)	
%     uLagr       = projectDG2LagrangeSub( uDG );
%     visualizeDataSub(g, uLagr, 'u', 'solu', k); 
    visualizeDataSub(g, particleTypeVector, 'particleType', 'solu', k); 
    end

    %bacLagr     = projectDG2LagrangeSub( uDGbac );
    %visualizeDataSub(g, bacLagr, 'bac', 'solBac', k);
    %biomVecLagr     = projectDG2LagrangeSub(bioConcVector);
    %visualizeDataSub(g, biomVecLagr, 'bioMass', 'bioMass', k);
    
%     particleDistribution = particleSizeDistribution(bulkVector.');
    


% printInfo(k,bulkVector,particle1List,particle2List,particle3List,particle4List,particle5List, particle6List, particle7List, particle8List, ...
%     particle9List, particle10List, particle11List, particle12List, particle13List, particle14List, particle15List, particle16List, particle17List, ...
%     particle18List, particle19List, particle20List)

%     numberOfSolidParticles = sum(particleDistribution(:,1) .* particleDistribution(:,2));

    if sum(bioMvector .* bulkVector ~= 0)
        error('NEEEEEEEEEEEIN!')
    end
%     particleDistribution
numEdgeTypes =  countEdgeTypes(g, bulkVector, POMVector, solidParticleList, ...
    edgeChargeVector, reactiveSurfaceVector, particleTypeVector);
POMsolidEdgeList = calculatePOMsolidEdgeList(g, bulkVector, POMVector, POMParticleList);
printInfoTUM(k,bulkVector,POMconcVector, concPOMAgent, edgeChargeVector, POMsolidEdgeList,...
    numEdgeTypes, totalPOMinputConc, sumExcessPOM, POMagentInput);
%     if(k < 25 || mod(k,25) == 0 || k == numOuterIt)
    if(k < 5 || mod(k,50) == 0 || k == numOuterIt)
        particleListHelper = particleList;
        particleList = solidParticleList;
        fileName    = ['FinalConfig/config','.', num2str(k),'.mat']; 
        save(fileName,'g','bulkVector','bulkTypeVector','POMconcVector', 'POMageVector', 'edgeChargeVector',...
            'POMVector','POMParticleList', 'particleList', 'reactiveSurfaceVector', 'particleTypeVector');
        particleList = particleListHelper;
    end
    
end  % for k
fclose( fileID );
% fclose( fileID_1);
end
