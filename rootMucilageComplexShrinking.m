function [rootGraph, bulkVector, rootVector, rootParticleList] = ...
    rootMucilageComplexShrinking(rootGraph, bulkVector, rootVector, rootParticleList,  rootDeadCellsInd)
    global notConnectedEdgesValue;
    for i = 1:numel(rootDeadCellsInd)
        
            neighInd = neighbors(rootGraph, rootDeadCellsInd(i));    
            rootNeighborCellsInd = neighInd(~ismember(neighInd,rootParticleList));
            
            h = rootDeadCellsInd(i) * ones(size(rootNeighborCellsInd,1),1);
            edgeInd = findedge(rootGraph, rootNeighborCellsInd, h);

            %nur van neumman weiterwachsen -> man dürfte nur Kanten mit 1 teilen
            if(rootGraph.Edges.Weight(edgeInd) < sqrt(2) * 1.1)
                rootGraph.Edges.Weight(edgeInd) = rootGraph.Edges.Weight(edgeInd).*notConnectedEdgesValue;
            end
            rootVector(rootDeadCellsInd(i)) = 0;
            bulkVector(rootDeadCellsInd(i)) = 0;
            rootParticleList = rootParticleList(rootParticleList ~= rootDeadCellsInd(i));
        
    end

end