function [bulkVector, MB_Vector, N_SVector, C_SVector, N_BVector, C_BVector ,C_MNVector, N_MNVector, MNVector, POMVector, C_POMconcVector, POMParticleList, POMageVector, CO2Vector,leakedNVector, C_EXTVector_t, C_PPlantVector,N_PPlantVector, C_PMNVector, N_PMNVector, CUE, ...
    f_C,f_BD_C,B_C,R, f_N,f_BD_N,B_N, R_leaked, f_POM_C, f_POM_N, f_MN_C, f_MN_N ] = ...
    calculate_C_N(g, parameters, bulkVector, MB_Vector, N_SVector, C_SVector, N_BVector, C_BVector, C_MNVector, N_MNVector, MNVector, POMVector, C_POMconcVector, reactiveSurfaceVector, POMParticleList, POMageVector,CO2Vector, leakedNVector,outerRootBorderInd, isMBfactor)
    
    dayinseconds = 24 * 60 * 60;% 
    numberoftstps = dayinseconds/parameters.tau_ode;
    range =5; % sqrt(1.94* 10^-8 cm^2/s*3600s)*10^4 (convert cm to microm) 
    C_EXTVector_t =  zeros(g.numT, 1);
    
   
    f_C = 0;
    f_BD_C = 0; 
    B_C = 0; 
    R = 0;

    f_N = 0;
    f_BD_N = 0; 
    B_N = 0; 
    R_leaked = 0;

    f_POM_C = 0;
    f_POM_N =0;
    f_MN_C = 0;
    f_MN_N = 0;
        
    RootExudates = 0;
    
    CO2Add = 0;
    %loop
    tolerance = 0.000001;
    for i = 1:numberoftstps%numberoftstps
        
    fprintf("step %d \n", i)
    leakedN_S = 0;

    C_PPlantVector_before = C_POMconcVector .* ~MNVector;
    N_PPlantVector_before = C_POMconcVector .* ~MNVector ./ parameters.C_N_POM;
    C_PMNVector_before = C_POMconcVector .* MNVector;
    N_PMNVector_before = C_POMconcVector .* MNVector / parameters.C_N_NM;

    
    previousConcentration_C =sum(C_MNVector+ C_SVector + C_BVector + CO2Vector+ C_POMconcVector);
    previousConcentration_N =sum(N_MNVector+ N_SVector + N_BVector + N_PMNVector_before + N_PPlantVector_before) + leakedN_S;

    before_C = sum(C_SVector + C_POMconcVector);     
   
    
    decaystep=tic;
    if(isMBfactor)
        MBfactor= 1 + sum(C_BVector)/50;
    else 
        MBfactor= 1;
    end
   MBfactor;
    [bulkVector,  POMVector, MNVector, C_POMconcVector, POMageVector, POMParticleList, C_SVector, N_SVector] = calculateonlyPOMdecay(g, parameters, bulkVector, POMVector,MNVector, C_POMconcVector,reactiveSurfaceVector, POMParticleList, POMageVector, C_SVector, N_SVector, MBfactor);
    fprintf('Time for POM MN Decay: %d \n', toc(decaystep))
    after_C = sum(C_SVector + C_POMconcVector);
    
    
    C_PPlantVector = C_POMconcVector .* ~MNVector;
    N_PPlantVector = C_POMconcVector .* ~MNVector ./ parameters.C_N_POM;
    C_PMNVector = C_POMconcVector .* MNVector;
    N_PMNVector = C_POMconcVector .* MNVector / parameters.C_N_NM;

    f_POM_C = f_POM_C + abs(sum(C_PPlantVector_before) - sum(C_PPlantVector));
    f_POM_N = f_POM_N + abs(sum(N_PPlantVector_before) - sum(N_PPlantVector));
    f_MN_C = f_MN_C + abs(sum(C_PMNVector_before) - sum(C_PMNVector));
    f_MN_N = f_MN_N + abs(sum(N_PMNVector_before) - sum(N_PMNVector));
        

    if(abs(before_C - after_C) > tolerance)
             abs(before_C - after_C)
              error('Falsch decay C %f', abs(before_C - after_C))
   end

    % Rootex_step = tic;
     sumC_S_before = sum(C_SVector);
     [C_SVector, N_SVector] = updateMucilage2(g, parameters,outerRootBorderInd, bulkVector,C_SVector, N_SVector);
     sumC_S = sum(C_SVector);
    % fprintf('Time for Rootex_step: %d \n', toc(Rootex_step))
    % sumC_S_avai = sum(C_SVector)
    RootExudates = RootExudates + (sumC_S - sumC_S_before);
    
     C_Nstep = tic;
     C_EXTVector =0;
     sumC_B_before = sum(C_BVector);
     sumC_S_before = sum(C_SVector);
     sumCO2_before = sum(CO2Vector);
     sumC_MN_before = sum(C_MNVector);

     sumN_B_before = sum(N_BVector);
     sumN_S_before = sum(N_SVector);
     sumN_MN_before = sum(N_MNVector);
     %TODO nur EPS produzieren wenn solid dran ist
  [ N_SVector, C_SVector, N_BVector, C_BVector,C_MNVector, N_MNVector, CO2Vector, C_EXTVector] = ...
    calculateMBSolid(g,parameters, N_SVector, C_SVector, N_BVector, C_BVector, C_MNVector, N_MNVector, CO2Vector);
     sumC_B = sum(C_BVector);
     sumC_S = sum(C_SVector);
     sumCO2 = sum(CO2Vector);
     sumC_MN = sum(C_MNVector);

     sumN_B = sum(N_BVector);
     sumN_S = sum(N_SVector);
     sumN_MN = sum(N_MNVector);
     
     f_C = f_C + abs(sumC_S - sumC_S_before);
     R = R + (sumCO2 - sumCO2_before);
     f_BD_C = f_BD_C + (sumC_MN - sumC_MN_before);
     B_C = B_C + (sumC_B - sumC_B_before);

     f_N = f_N + abs(sumN_S - sumN_S_before);
    % hallo = f_C/f_N
     f_BD_N = f_BD_N + (sumN_MN - sumN_MN_before);
     B_N = B_N + (sumN_B - sumN_B_before);

     if(isnan(sumC_S))
         falsch = 1;
     end
     
     CO2Add = CO2Add + abs(sumCO2_before - sumCO2);
     CUE = (f_C - CO2Add)/f_C
     
    % if(numel(find(C_SVector < 0)) > 0)
    %     error('Falsch C_S')
    % end
    sumC_B_after = sum(C_BVector);
    fprintf('Time for MB step: %d \n', toc(C_Nstep))
    C_EXTVector_t = C_EXTVector_t + C_EXTVector;
    
    %% Diff and spread
    spreadstep=tic;
    maxValues = ones(g.numT, 1).* parameters.maxConcC_B - C_MNVector * parameters.waterContentNecromass;
    maxValues(maxValues < 0) = 0;
    C_BVector = spreadConcentration(g, C_BVector, bulkVector, maxValues, ones(g.numT, 1).*parameters.minConC_B);
    N_BVector = C_BVector ./ parameters.C_N_B;
    fprintf('Time for spreadstp: %d \n', toc(spreadstep))



    f_BD_C = f_BD_C + sum(C_BVector .* (C_BVector < parameters.minConC_B));
    f_BD_N = f_BD_N + sum(N_BVector .* (C_BVector < parameters.minConC_B));


    C_SVector = C_SVector + C_BVector .* (C_BVector < parameters.minConC_B); 
    C_BVector(C_BVector < parameters.minConC_B) = 0;

    N_SVector = N_SVector + N_BVector .* (C_BVector < parameters.minConC_B); 
    N_BVector(C_BVector < parameters.minConC_B) = 0;
    
    MB_Vector =  C_BVector >= parameters.minConC_B;
    n_MB = sum(MB_Vector)
    %MN not part of MB -> POM
   % C_BVector(:) = 0;
    indx = intersect(find(C_BVector == 0), find(C_MNVector > 0))
    len_POM_List  = length(POMParticleList)
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
        %len_POM_List  = length(POMParticleList)
       %totalPOMinputConc = totalPOMinputConc + length(globalIndNewPOMparticle);
    end
   % sumMNVector = sum(MNVector);

    


    diffstep = tic;
    %sumC_S_before = sum(C_SVector)
    C_SVector = easyDiffusiveStep(g, C_SVector, bulkVector, MB_Vector, range);
    %sumC_S_after= sum(C_SVector)
    N_SVector = easyDiffusiveStep(g, N_SVector, bulkVector, MB_Vector, range);
    leakedN_S = sum(N_SVector .* parameters.N_leakage);
    R_leaked= R_leaked  + leakedN_S;

    N_SVector = N_SVector .* (1-parameters.N_leakage);
    leakedNVector = leakedNVector + N_SVector .* parameters.N_leakage;
    %N_SVector = C_SVector ./ parameters.C_N_S;
    fprintf('Time for diffstep: %d \n', toc(diffstep))

    currentConcentration_C = sum(C_MNVector+ C_SVector + C_BVector + CO2Vector +C_POMconcVector);
    % fprintf("abs(previousConcentration_C - currentConcentration_C) %f \n", abs(previousConcentration_C - currentConcentration_C))
   if(abs(previousConcentration_C - currentConcentration_C) > tolerance)
             %abs(previousConcentration_C - currentConcentration_C)
              fprintf('Falsch C:N C %f \n', abs(previousConcentration_C - currentConcentration_C))
   end

    C_PPlantVector = C_POMconcVector .* ~MNVector;
    N_PPlantVector = C_POMconcVector .* ~MNVector ./ parameters.C_N_POM;
    C_PMNVector = C_POMconcVector .* MNVector;
    N_PMNVector = C_POMconcVector .* MNVector / parameters.C_N_NM;



   currentConcentration_N = sum(N_MNVector+ N_SVector + N_BVector + N_PMNVector + N_PPlantVector) + leakedN_S;
   %fprintf("abs(previousConcentration_N - currentConcentration_N) %f \n", abs(previousConcentration_N - currentConcentration_N))
   if(abs(previousConcentration_N - currentConcentration_N) > tolerance)
            % abs(previousConcentration_N - currentConcentration_N)
             fprintf('Falsch C:N  N %f \n', abs(previousConcentration_N - currentConcentration_N))

   end


    %Concentration_Phase = find(C_SVector < 0)
%numel(find(C_BVector > 0))
%C_BVector(find(C_BVector > 0))
    
    end % end loop
    CUE = (f_C - CO2Add)/f_C
end