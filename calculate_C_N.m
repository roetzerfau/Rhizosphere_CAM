function [bulkVector, MB_Vector, N_SVector, C_SVector, N_BVector, C_BVector ,C_MNVector, N_MNVector, MNVector, POMVector, C_POMconcVector, POMParticleList, POMageVector, C_EXTVector_t] = ...
    calculate_C_N(g, parameters, bulkVector, MB_Vector, N_SVector, C_SVector, N_BVector, C_BVector, C_MNVector, N_MNVector, MNVector, POMVector, C_POMconcVector, reactiveSurfaceVector, POMParticleList, POMageVector,outerRootBorderInd)

    dayinseconds = 24 * 60 * 60;
    numberoftstps = dayinseconds/parameters.tau_ode;
    range =round( 5/(60*60) * parameters.tau_ode);
    C_EXTVector_t =  zeros(g.numT, 1);
    for i = 1:numberoftstps
        fprintf("step %d \n", i)
    %before = sum(C_SVector)     
    [bulkVector,  POMVector, MNVector, C_POMconcVector, POMageVector, POMParticleList, C_SVector, N_SVector] = calculateonlyPOMdecay(g, parameters, bulkVector, POMVector,MNVector, C_POMconcVector,reactiveSurfaceVector, POMParticleList, POMageVector, C_SVector, N_SVector);
    %after = sum(C_SVector)
    
    Rootex_step = tic;
    [C_SVector, N_SVector] = updateMucilage2(g, parameters,outerRootBorderInd, bulkVector,C_SVector, N_SVector);
    fprintf('Time for Rootex_step: %d \n', toc(Rootex_step))
    %sumC_S_avai = sum(C_SVector)
    
     C_Nstep = tic;
     sumC_B_before = sum(C_BVector)
    [ N_SVector, C_SVector, N_BVector, C_BVector,C_MNVector, N_MNVector, C_EXTVector] = ...
    calculateMBSolid(g,parameters, N_SVector, C_SVector, N_BVector, C_BVector, C_MNVector, N_MNVector);
    %sumC_S = sum(C_SVector)
    sumC_B_after = sum(C_BVector)
    fprintf('Time for C_Nstep: %d \n', toc(C_Nstep))
    C_EXTVector_t = C_EXTVector_t + C_EXTVector;
    
    %% Diff and spread
    spreadstep=tic;
    maxValues = ones(g.numT, 1).* parameters.maxConcC_B - C_MNVector * parameters.waterContentNecromass;
    maxValues(maxValues < 0) = 0;
    C_BVector = spreadConcentration(g, C_BVector, bulkVector, maxValues, ones(g.numT, 1).*parameters.minConC_B);
    N_BVector = C_BVector ./ parameters.C_N_B;
    fprintf('Time for spreadstp: %d \n', toc(spreadstep))



    
    C_SVector = C_SVector + C_BVector .* (C_BVector < parameters.minConC_B); 
    C_BVector(C_BVector < parameters.minConC_B) = 0;
    
    N_SVector = N_SVector + N_BVector .* (C_BVector < parameters.minConC_B); 
    N_BVector(C_BVector < parameters.minConC_B) = 0;
    
    MB_Vector =  C_BVector >= parameters.minConC_B;

    %MN not part of MB -> POM
    indx = intersect(find(C_BVector == 0), find(C_MNVector > 0));
    for i = 1:numel(indx)
        globalIndNewPOMparticle = indx(i);

        MNVector(globalIndNewPOMparticle) = 1;
        if(bulkVector(globalIndNewPOMparticle) == 1 || POMVector(globalIndNewPOMparticle) == 1)
            fprintf('falsch')
        end
        bulkVector(globalIndNewPOMparticle) = 1;
        POMVector(globalIndNewPOMparticle) = 1;
        C_POMconcVector(globalIndNewPOMparticle) = C_POMconcVector(globalIndNewPOMparticle) + C_MNVector(globalIndNewPOMparticle);
        C_MNVector(globalIndNewPOMparticle) = 0;
        N_MNVector(globalIndNewPOMparticle) = 0;
        POMageVector(globalIndNewPOMparticle) = 1;
        POMParticleList{length(POMParticleList) + 1} = globalIndNewPOMparticle;
       %totalPOMinputConc = totalPOMinputConc + length(globalIndNewPOMparticle);
    end
    sumMNVector = sum(MNVector)

    


    diffstep = tic;
    %sumC_S_before = sum(C_SVector)
    C_SVector = easyDiffusiveStep(g, C_SVector, bulkVector, range);
    %sumC_S_after= sum(C_SVector)
    N_SVector = easyDiffusiveStep(g, N_SVector, bulkVector, range);
    N_SVector = N_SVector .* (1-parameters.N_leakage);
    %N_SVector = C_SVector ./ parameters.C_N_S;
    fprintf('Time for diffstep: %d \n', toc(diffstep))

%numel(find(C_BVector > 0))
%C_BVector(find(C_BVector > 0))
    end
end