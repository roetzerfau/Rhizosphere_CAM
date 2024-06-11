 function [bulkVector,  POMVector, C_POMconcVector, POMageVector, POMParticleList, C_SVector, N_SVector] = calculateonlyPOMdecay(g, parameters, bulkVector, POMVector, C_POMconcVector,reactiveSurfaceVector, POMParticleList, POMageVector, C_SVector, N_SVector)
    
 numFluidNeighVector = calculateNumFluidNeighbors(g, bulkVector, reactiveSurfaceVector, POMVector, 2);
    POMsolidEdgeList = calculatePOMsolidEdgeList(g, bulkVector, POMVector, POMParticleList);
    
   
    % update POM conc
    POMconcVectorOld = C_POMconcVector;
    dayinsecondsec = 24 * 60 * 60;
    tau = 1/dayinsecondsec * parameters.tau_ode;
    
    %tau = parameters.tau_ode;
    %parameters.POMdecayRate = 2.89 * 10^-6;
    % only POM particles that are attached to reactive solid surface are
    % decaying
    for i = 1 : length(POMParticleList)
        if (sum(reactiveSurfaceVector(POMsolidEdgeList{i})) > 0)... % + sum(edgeChargeVector(POMsolidEdgeList{i}))
                && (sum(numFluidNeighVector(POMParticleList{i})) > 0)
            concOld = sum(C_POMconcVector(POMParticleList{i}));
            particleDecayRate = parameters.POMdecayRate * sum(numFluidNeighVector(POMParticleList{i})) / ...
                (sum(numFluidNeighVector(POMParticleList{i})) + length(POMsolidEdgeList{i}));
            concNew = concOld * exp(- particleDecayRate * tau);
            concDiff = concOld - concNew;
            C_POMconcVector(POMParticleList{i}) = C_POMconcVector(POMParticleList{i}) - concDiff * ...
                numFluidNeighVector(POMParticleList{i}) / sum(numFluidNeighVector(POMParticleList{i}));
            
            C_SVector(POMParticleList{i}) = C_SVector(POMParticleList{i}) +concDiff;
            N_SVector(POMParticleList{i}) = (N_SVector(POMParticleList{i}) + concDiff) / parameters.C_N_POM; 
        
        end
    end
    
    % set POM conc. below threshold to zero
    C_POMconcVector(C_POMconcVector < parameters.POMminConc) = 0;
    
    POMdecayVector = POMconcVectorOld - C_POMconcVector;
    

    % save amount of decayed POM
    decayedPOMfromParticle = zeros(length(POMParticleList),1);
    for i = 1 : length(POMParticleList)
        decayedPOMfromParticle(i) = sum(POMdecayVector(POMParticleList{i}));
    end
    
    
        % update bulkVector etc.
    indPOMchanged = find( ( POMVector == 1 ) &  ( C_POMconcVector == 0 ) );
    
    POMVector(indPOMchanged) = 0;
    bulkVector(indPOMchanged) = 0;
    POMageVector(indPOMchanged) = 0;
    
    % Update POM age
    POMageVector(POMageVector > 0) = POMageVector(POMageVector > 0) + 1;

    
    % if POM cells dis                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             appeared
    if (~isempty(indPOMchanged))
%         [POMParticleList, ~] = createPOMParticleList(POMVector);
        % remove degraded POM cells from POMParticleList
        for i = 1 : length(POMParticleList)
            POMParticleList{i}(find(ismember(POMParticleList{i}, indPOMchanged))) = [];
        end
        % remove all POM particles that are completely degraded from POMParticleList
        for i = length(POMParticleList) : -1 : 1
            if (length(POMParticleList{i})) == 0
                POMParticleList(i) = [];
            end
        end
        
        % find POM particles that are separated into multiple parts due to
        % degradation
        POMParticleListOld = POMParticleList;
        POMparticlesToRemove = [];
        for i = 1 : length(POMParticleListOld)
            helperVector = zeros(size(POMVector));
            helperVector(POMParticleListOld{i}) = 1;
            % if particle was separated, length(POMListHelper) > 1
            % check periodicity...
            [POMListHelper, ~] = createPOMParticleListNew(helperVector);
            if length(POMListHelper) > 1
                % add separated POM particles at end and remove original
                % from POMParticleList
                POMParticleList = [POMParticleList POMListHelper];
                POMparticlesToRemove = [POMparticlesToRemove i];
            end
        end
        POMParticleList(POMparticlesToRemove) = [];
        
    end
    
    
 end
 
