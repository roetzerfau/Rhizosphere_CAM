function [mucilageVector, mucilageEPSConcVector, mucilageSurfaceVector, MucilageagentAge] = updateMucilageEPS(g,parameters,outerRootBorderInd, EPSInd, bulkVector,mucilageVector, mucilageEPSConcVector, MucilageagentAge, mucilageSurfaceVector)
    

    if(parameters.mucilageGrowing == 1)
        nA = sum(bulkVector(outerRootBorderInd) == 0);
        constA = p.constantMucilageDeposition/nA;
            for i = 1:numel(outerRootBorderInd)
               if(bulkVector(outerRootBorderInd(i)) == 0)
                mucilageEPSConcVector(outerRootBorderInd(i)) = mucilageEPSConcVector(outerRootBorderInd(i)) + constA;           
               end
            end
    end 
    
    
    
    
    
    
    mucilageEPSConcVector = spreadConcentration(g, mucilageEPSConcVector, bulkVector, p.normalMucilageConcentration, p.minConcMucilage);
    
    mucilageVector = mucilageEPSConcVector > p.minConcMucilage;
    
    
    
    mucilageParticleList = cell(1,1);
    mucilageParticleList{1} = find(mucilageVector == 1);
    mucilageSolidEdgeList = calculatePOMsolidEdgeList(g, bulkVector, mucilageVector, mucilageParticleList);
    MucilageagentAge(mucilageSolidEdgeList{1}) = 0;
    mucilageSurfaceVector(mucilageSolidEdgeList{1}) = 1;

    mucilageTest = (0*ones( g.numCE , 1 ));   
    allrellevantEdges = g.CE0T(find(bulkVector == 1),:);
    mucilageTest(allrellevantEdges) = 1;
    mucilageSurfaceVector = mucilageSurfaceVector .* mucilageTest;
    

end