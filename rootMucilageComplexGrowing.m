function [rootComplexGraph, bulkVector, rootComplexList, rootPressureDistributionVector] = ...
rootMucilageComplexGrowing(g, rootComplexGraph, bulkVector, rootComplexVector,  rootComplexList, amountNewCells)
    
    N = g.NX;
    rootSourceCell = rootComplexList(1);
    notConnectedEdgesValue = N* N * 2;
   
    %% Calculate possible Cells for Growing
    
    freeCellsInd = find((rootComplexVector) ~= 1);
    
    [TR,d] = shortestpathtree(rootComplexGraph,freeCellsInd,rootSourceCell);
    
    [sortedd, I] = sort(d);
    
    freeCellsInd = freeCellsInd(I);
    
    %
    i = 1;
    newCellsInd = [];
    while(amountNewCells > numel(newCellsInd) && numel(freeCellsInd) >= i)
            st = stencil(g.NX,g.NX,freeCellsInd(i),1);
            %look if continuous growing without holes or skipped
            %connections
            isCellNeighboredRoot = sum(ismember(st(2:end),rootComplexList),'all') > 0;
            if(isCellNeighboredRoot == 0)               
                %fprintf('geht nicht')
            end

            if(bulkVector(freeCellsInd(i)) == 0 && isCellNeighboredRoot)
                neighInd = neighbors(rootComplexGraph, freeCellsInd(i));    
                rootBorderCellsInd = intersect(neighInd,rootComplexList);

                h = freeCellsInd(i) * ones(size(rootBorderCellsInd,1),1);
                edgeInd = findedge(rootComplexGraph, rootBorderCellsInd, h);

                %nur van neumman weiterwachsen -> man dürfte nur Kanten mit 1 teilen
                %TODO sicher stellen dass überall richtige Werte
                if(rootComplexGraph.Edges.Weight(edgeInd) >= sqrt(2) * 1.1)
                    rootComplexGraph.Edges.Weight(edgeInd) = rootComplexGraph.Edges.Weight(edgeInd)./notConnectedEdgesValue;
                end
                bulkVector(freeCellsInd(i)) = 1;
                rootComplexList = [rootComplexList freeCellsInd(i)];
                newCellsInd = [newCellsInd freeCellsInd(i)];
            else
                %fprintf('Cell cant grow \n');
            end
            i = i+1;      
    end 

    
    
    
    %% Cells which couldnt grow are converted to pressure points
    visitedFreeCells = 1:i-1;
    if(numel(newCellsInd) == 0)
        pressurePointsInd = freeCellsInd;
    else
        growedFreeCells = sum(freeCellsInd == newCellsInd,2);
        %pressurePointsInd = freeCellsInd(visitedFreeCells ~= growedFreeCells);
        pressurePointsInd = freeCellsInd(visitedFreeCells);
        pressurePointsInd = freeCellsInd(~growedFreeCells(visitedFreeCells));
    end
    pressurePointsConnectedBulk = 0 * ones(g.numT, 1);
    rootPressureDistributionVector = 0 * ones(g.numT, 1);
    %die zellen die nicht wachsen konnten druckpunkte
    %die freien cellen aus array rausmachen, die nächsten als Druckpunkte
    %Mit pressurePoints verbunden
    bulkVector_bw = (reshape(bulkVector, [N N]));
    CC = bwconncomp(bulkVector_bw, 4);
    pressurePointsConnectedBulk(:)= 0;
    pressurePointsConnectBulkInd = [];
    for p = 1:numel(pressurePointsInd)
        oo = cellfun(@(m) m == pressurePointsInd(p),CC.PixelIdxList, 'UniformOutput', false);
        Match1 = cellfun(@sum, oo);
        Match = find(Match1 == 1);
        
        if(~isempty(Match))
            pressurePointsConnectedBulkInd  =  CC.PixelIdxList{Match};    
            pressurePointsConnectedBulk(pressurePointsConnectedBulkInd) = 1;
            pressurePointsConnectBulkInd = [pressurePointsConnectBulkInd pressurePointsInd];
        end
    end
     %-------------------------------------------------------
    % Pressure Distributoin - pressure Points ohne Verbindung keine
    % Auswirkung
    [X,Y] = meshgrid(1:N, 1:N);
    rootPressureDistributionVector(:) = 0;
    rootPressureDistributionVector = reshape(rootPressureDistributionVector, [N N]);
    for i = 1:numel(pressurePointsConnectBulkInd)
        
        
        X_0 = X(pressurePointsConnectBulkInd(i));
        Y_0 = Y(pressurePointsConnectBulkInd(i));
        F = -((X-X_0).^2 + (Y-Y_0).^2) + (N/4)^2 ;
        F(F < 0) = 0;
        rootPressureDistributionVector = rootPressureDistributionVector + F;
    end
    
    rootPressureDistributionVector = reshape(rootPressureDistributionVector, [N*N 1]);
    rootPressureDistributionVector = normalize(rootPressureDistributionVector,'range',[0 1]);
    %---------------------------------------------------------------------------------
    % Combine Pressure mit Connected with Root
    pressureArea = pressurePointsConnectedBulk;
    pressureArea(pressurePointsInd) = 1;
    pressureArea(rootComplexList) = 0;
    rootPressureDistributionVector = rootPressureDistributionVector.*pressureArea;
    
    
    %todo
    %rootPressureDistributionVector(rootPressureDistributionVector > 0) = normalize(rootPressureDistributionVector(rootPressureDistributionVector > 0),'range',[0 1]);
    
    
    
    
    
    




%pressureCells werden zu potentiellen neuen Zellen
% da wo rootPressuredistribution == 1 ,da neue Zellen
end