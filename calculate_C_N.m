function [bulkVector, MB_Vector, N_SVector, C_SVector, N_BVector, C_BVector ,C_MNVector, N_MNVector, MNVector, POMVector, C_POMconcVector, POMParticleList, POMageVector, CO2Vector, C_EXTVector_t, C_PPlantVector,N_PPlantVector, C_PMNVector, N_PMNVector, sumleakedN_S, CUE ] = ...
    calculate_C_N(g, parameters, bulkVector, MB_Vector, N_SVector, C_SVector, N_BVector, C_BVector, C_MNVector, N_MNVector, MNVector, POMVector, C_POMconcVector, reactiveSurfaceVector, POMParticleList, POMageVector,CO2Vector, outerRootBorderInd)
    
    dayinseconds = 24 * 60 * 60;% 
    numberoftstps = dayinseconds/parameters.tau_ode;
    range =round( 5/(60*60) * parameters.tau_ode);
    C_EXTVector_t =  zeros(g.numT, 1);
    sumleakedN_S = 0;

    for i = 1:numberoftstps
    fprintf("step %d \n", i)
    leakedN_S = 0;

    C_PPlantVector = C_POMconcVector .* ~MNVector;
    N_PPlantVector = C_POMconcVector .* ~MNVector ./ parameters.C_N_POM;
    C_PMNVector = C_POMconcVector .* MNVector;
    N_PMNVector = C_POMconcVector .* MNVector / parameters.C_N_NM;

    
    previousConcentration_C =sum(C_MNVector+ C_SVector + C_BVector + CO2Vector+ C_POMconcVector);
    previousConcentration_N =sum(N_MNVector+ N_SVector + N_BVector + N_PMNVector + N_PPlantVector) + leakedN_S;

    %before = sum(C_SVector + C_POMconcVector)     
    [bulkVector,  POMVector, MNVector, C_POMconcVector, POMageVector, POMParticleList, C_SVector, N_SVector] = calculateonlyPOMdecay(g, parameters, bulkVector, POMVector,MNVector, C_POMconcVector,reactiveSurfaceVector, POMParticleList, POMageVector, C_SVector, N_SVector);
    %after = sum(C_SVector + C_POMconcVector)
    
    %Rootex_step = tic;
    %[C_SVector, N_SVector] = updateMucilage2(g, parameters,outerRootBorderInd, bulkVector,C_SVector, N_SVector);
    %fprintf('Time for Rootex_step: %d \n', toc(Rootex_step))
    %sumC_S_avai = sum(C_SVector)
    
     C_Nstep = tic;
     C_EXTVector =0;
     sumC_B_before = sum(C_BVector);
     sumC_S_before = sum(C_SVector);
     sumCO2_before = sum(CO2Vector);
     %TODO nur EPS produzieren wenn solid dran ist
  [ N_SVector, C_SVector, N_BVector, C_BVector,C_MNVector, N_MNVector, CO2Vector, C_EXTVector] = ...
    calculateMBSolid(g,parameters, N_SVector, C_SVector, N_BVector, C_BVector, C_MNVector, N_MNVector, CO2Vector);
     sumC_S = sum(C_SVector);
     if(isnan(sumC_S))
         falsch = 1;
     end
     U = abs(sumC_S_before - sumC_S);
     sumCO2 = sum(CO2Vector);
     CO2Add = abs(sumCO2_before - sumCO2);
     CUE = (U - CO2Add)/U;
     
    if(numel(find(C_SVector < 0)) > 0)
        error('Falsch C_S')
    end
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
    leakedN_S = sum(N_SVector .* parameters.N_leakage);
    sumleakedN_S= sumleakedN_S  + leakedN_S;
    N_SVector = N_SVector .* (1-parameters.N_leakage);
   
    %N_SVector = C_SVector ./ parameters.C_N_S;
    fprintf('Time for diffstep: %d \n', toc(diffstep))

    currentConcentration_C = sum(C_MNVector+ C_SVector + C_BVector + CO2Vector +C_POMconcVector);
   if(abs(previousConcentration_C - currentConcentration_C) > 0.1)
             abs(previousConcentration_C - currentConcentration_C)
              error('Falsch C:N C %f', abs(previousConcentration_C - currentConcentration_C))
   end

    C_PPlantVector = C_POMconcVector .* ~MNVector;
    N_PPlantVector = C_POMconcVector .* ~MNVector ./ parameters.C_N_POM;
    C_PMNVector = C_POMconcVector .* MNVector;
    N_PMNVector = C_POMconcVector .* MNVector / parameters.C_N_NM;



   currentConcentration_N = sum(N_MNVector+ N_SVector + N_BVector + N_PMNVector + N_PPlantVector) + leakedN_S;
   if(abs(previousConcentration_N - currentConcentration_N) > 0.1)
             abs(previousConcentration_N - currentConcentration_N)
              error('Falsch C:N  N %f', abs(previousConcentration_N - currentConcentration_N))

   end


    %Concentration_Phase = find(C_SVector < 0)
%numel(find(C_BVector > 0))
%C_BVector(find(C_BVector > 0))
    
    end
end