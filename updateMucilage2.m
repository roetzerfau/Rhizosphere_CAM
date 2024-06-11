function [C_SVector, N_SVector] = updateMucilage2(g, parameters,outerRootBorderInd, bulkVector,C_SVector, N_SVector)


    if(parameters.mucilageGrowing == 1)
        %nA = sum(bulkVector(outerRootBorderInd) == 0);
        %constA = p.constantMucilageDeposition/nA;
            for i = 1:numel(outerRootBorderInd)
               if(bulkVector(outerRootBorderInd(i)) == 0)
                %mucilageConcVector(outerRootBorderInd(i)) = mucilageConcVector(outerRootBorderInd(i)) + constA;
                C_SVector(outerRootBorderInd(i)) = C_SVector(outerRootBorderInd(i)) + parameters.startConcPOM/50;
                N_SVector(outerRootBorderInd(i)) = N_SVector(outerRootBorderInd(i)) /parameters.C_N_Root;
               end
            end
    end 

end