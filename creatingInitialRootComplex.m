function [rootVector, mucilageVector, mucilageConcVector, mucilageSurfaceVector, rootComplexList, rootComplexGraph, mucilageGraph] = ...
creatingInitialRootComplex(g, bulkVector)

    N = g.NX;
    notConnectedEdgesValue = N*N * 2;
    rootVector = 0 * ones(g.numT, 1);
    rootComplexList = [];
    mucilageVector = 0 * ones(g.numT, 1);
    mucilageConcVector = 0 * ones(g.numT, 1);
    mucilageSurfaceVector = (0*ones( g.numCE , 1 ));   

    diagVec1 = sparse(repmat([ones(N-1, 1); 0], N, 1));  % Make the first diagonal vector
                                                 %   (for horizontal connections)
    diagVec1 = diagVec1(1:end-1);                % Remove the last value
    diagVec2 = sparse([0; diagVec1(1:(N*(N-1)))] * sqrt(2));       % Make the second diagonal vector
                                                 %   (for anti-diagonal connections)
    diagVec3 = sparse(ones(N*(N-1), 1));                 % Make the third diagonal vector
                                                 %   (for vertical connections)
    diagVec4 = sparse(diagVec2(2:end-1));                % Make the fourth diagonal vector
                                                 %   (for diagonal connections)
    adj = diag(diagVec1, 1)+...                  % Add the diagonals to a zero matrix
          diag(diagVec2, N-1)+...
          diag(diagVec3, N)+...
          diag(diagVec4, N+1);
    adj = adj+adj.';                             % Add the matrix to a transposed copy of
                                                 %   itself to make it symmetric

    adj = adj .*notConnectedEdgesValue;

    rootComplexGraph = graph(adj);
    mucilageGraph = rootComplexGraph;
    clear adj diagVec1 diagVec2 diagVec3 diagVec4

    %find nearest pore space to center of domain 
    centerOfDomain = g.NX/2;
    rootIntialCellInd = centerOfDomain * g.NX + centerOfDomain;
    
    [TR,d] = shortestpathtree(rootComplexGraph,rootIntialCellInd);
    [sortedd, I] = sort(d);
    j =find(bulkVector(I) == 0,1);
    rootNewCellsInd = I(j);
    rootVector(rootNewCellsInd) = 1;
    rootComplexList = [rootComplexList rootNewCellsInd];
    
end