function [rootComplexGraph, bulkVector, rootVector, mucilageVector, rootComplexParticleList, rootPressureDistributionVector] = ...
    rootMucilageComplexShrinking(g, parameters, rootComplexGraph, bulkVector, rootVector, mucilageVector, rootComplexParticleList,  rootSourceCell, rootNr)
    

 [mucilageConcVector, rootVector, rootParticleList] = calculateMucilageDecay(g, parameters, bulkVector, rootVector, rootParticleList, mucilageConcVector);
    rootDeadCellsInd = find( ( rootVector == 1 ) &  ( mucilageConcVector == 0 ) );
    
    for i = 1:numel(rootDeadCellsInd)
        
            neighInd = neighbors(rootGraph, rootDeadCellsInd(i));    
            rootNeighborCellsInd = neighInd(~ismember(neighInd,rootParticleList{rootNr}));
            
            h = rootDeadCellsInd(i) * ones(size(rootNeighborCellsInd,1),1);
            edgeInd = findedge(rootGraph, rootNeighborCellsInd, h);

            %nur van neumman weiterwachsen -> man dürfte nur Kanten mit 1 teilen
            if(rootGraph.Edges.Weight(edgeInd) < sqrt(2) * 1.1)
                rootGraph.Edges.Weight(edgeInd) = rootGraph.Edges.Weight(edgeInd).*notConnectedEdgesValue;
            end
            rootVector(rootDeadCellsInd(i)) = 0;
            bulkVector(rootDeadCellsInd(i)) = 0;
            rootParticleList{rootNr} = rootParticleList{rootNr}(rootParticleList{rootNr} ~= rootDeadCellsInd(i));
        
    end

end