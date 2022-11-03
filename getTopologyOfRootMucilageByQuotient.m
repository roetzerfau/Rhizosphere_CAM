function [rootVector mucilageVector] = getTopologyOfRootMucilageByQuotient(g,rootComplexGraph,rootComplexList, q_mucilage)
rootVector = 0 * ones(g.numT, 1);
mucilageVector = 0 * ones(g.numT, 1);
    
rootMucilageComplexVector = zeros(g.numT, 1);
rootMucilageComplexVector(rootComplexList) = 1;
% rootVector = rootMucilageComplexVector;
[relMap, borderpointsVector] = calculateRelativeDistanceMap(g,rootComplexGraph,rootComplexList,rootMucilageComplexVector ,rootComplexList(1));
relMap(relMap < 0) = -1;%fürs abspeichern
a = sum(relMap >= 0) == numel(rootComplexList);
% das ist manchmal das problem, könnte an periodizutät liegen


[relMap_sorted, I] = sort(relMap, 'descend');
amountMucilageCells = ceil(numel(rootComplexList) * q_mucilage);
mucilageVector(I(1:amountMucilageCells)) = 1;
rootVector(I(amountMucilageCells+1:numel(rootComplexList))) = 1;
end