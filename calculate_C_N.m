function [MB_Vector, N_SVector, C_SVector, N_BVector, C_BVector , MNVector, POMVector, C_POMconcVector, POMParticleList, POMageVector] = ...
    calculate_C_N(g, parameters, bulkVector, MB_Vector, N_SVector, C_SVector, N_BVector, C_BVector , MNVector, POMVector, C_POMconcVector, reactiveSurfaceVector, POMParticleList, POMageVector,outerRootBorderInd)

    dayinseconds = 24 * 60 * 60;
    numberoftstps = dayinseconds/parameters.tau_ode;
    range =round( 5/(60*60) * parameters.tau_ode);
    
    for i = 1:numberoftstps
        
    before = sum(C_SVector)     
    [bulkVector,  POMVector, C_POMconcVector, POMageVector, POMParticleList, C_SVector, N_SVector] = calculateonlyPOMdecay(g, parameters, bulkVector, POMVector, C_POMconcVector,reactiveSurfaceVector, POMParticleList, POMageVector, C_SVector, N_SVector);
    after = sum(C_SVector)
%     
    Rootex_step = tic;
    [C_SVector, N_SVector] = updateMucilage2(g, parameters,outerRootBorderInd, bulkVector,C_SVector, N_SVector);
    fprintf('Time for Rootex_step: %d \n', toc(Rootex_step))
    sumC_S_avai = sum(C_SVector)
    
     C_Nstep = tic;
    [ N_SVector, C_SVector, N_BVector, C_BVector,MNVector] = ...
    calculateMBSolid(g,parameters, N_SVector, C_SVector, N_BVector, C_BVector, bulkVector,  MNVector);
    sumC_S = sum(C_SVector)
    sumC_B = sum(C_BVector)
    fprintf('Time for C_Nstep: %d \n', toc(C_Nstep))

    
    %% Diff and spread
    spreadstep=tic;
    C_BVector = spreadConcentration(g, C_BVector, bulkVector,parameters.maxConcC_B, parameters.minConC_B);
    N_BVector = C_BVector ./ parameters.C_N_B;
    fprintf('Time for spreadstp: %d \n', toc(spreadstep))
    
    C_SVector = C_SVector + C_BVector .* (C_BVector < parameters.minConC_B); 
    C_BVector(C_BVector < parameters.minConC_B) = 0;
    
    N_SVector = N_SVector + N_BVector .* (C_BVector < parameters.minConC_B); 
    N_BVector(C_BVector < parameters.minConC_B) = 0;
    
    MB_Vector =  C_BVector > parameters.minConC_B;
    
    diffstep = tic;
    C_SVector = easyDiffusiveStep(g, C_SVector, bulkVector, range);
    N_SVector = easyDiffusiveStep(g, N_SVector, bulkVector, range);
    %N_SVector = C_SVector ./ parameters.C_N_S;
    fprintf('Time for diffstep: %d \n', toc(diffstep))

%numel(find(C_BVector > 0))
%C_BVector(find(C_BVector > 0))
    end
end