function [rootComplexGraph, bulkVector, rootComplexList, rootPressureDistributionVector , optimalFutureRootComplexList] = ...
rootMucilageComplexGrowing(g, rootComplexGraph, bulkVector, rootComplexVector,  rootComplexList, amountNewCells)
    
    N = g.NX;
    rootSourceCell = rootComplexList(1);
    notConnectedEdgesValue = N* N * 2;
    newCellsInd = [];
    %% Calculate possible Cells for Growing
    % übereinstimung amountCell von freeCells mit borderpoints 
    %diese wachsen, rest Druckpunkte
    %neue Borderpunkte
    % flaggen  connect mit Borderpoint pore -> wachsen
    %          connect mit Borderpoint bulk -> druckpunkt   
    %nach jedem durchgang updaten
    %
    %
    %Vielleicht nach jedem Wachstumsschritt move Particle von Partivle
    %direkt am Wachstum dran 
    for i = 1:amountNewCells
        % Suchen
        freeCellsInd = find((rootComplexVector) ~= 1);
        
        [TR,d] = shortestpathtree(rootComplexGraph,freeCellsInd,rootSourceCell);

        [sortedd, I] = sort(d);
        freeCellsInd = freeCellsInd(I);
        
        outerborder = sortedd < notConnectedEdgesValue * 1.1;%so könnte man theoretisch auch weit entferntest border pixel findne
        outerborderInd = freeCellsInd(outerborder);
        %outerborderInd(1) wenn da particle dran, dann moveParticle
        %pro Zeitschritt:
        % - es soll nur an der Wurzelbroder was wachsen
        % es soll bei jedem wachsen, etwas weggestoßen werdne wenn nötig
        occupiedCellsInd = outerborderInd(bulkVector(outerborderInd) == 1);
        
        
        
        newCellInd = outerborderInd(bulkVector(outerborderInd) == 0);
        if(numel(newCellInd) ~= 0)
            newCellInd = newCellInd(1);
        else
            %newCellInd = freeCellsInd(end);
            newCellInd = [];
            break;
        end
        
        % Markieren
        neighInd = neighbors(rootComplexGraph, newCellInd);    
        rootBorderCellsInd = intersect(neighInd,rootComplexList);

        h = newCellInd * ones(size(rootBorderCellsInd,1),1);
        edgeInd = findedge(rootComplexGraph, rootBorderCellsInd, h);

        %nur van neumman weiterwachsen -> man dürfte nur Kanten mit 1 teilen
        %TODO sicher stellen dass überall richtige Werte
        if(rootComplexGraph.Edges.Weight(edgeInd) >= sqrt(2) * 1.1)
            rootComplexGraph.Edges.Weight(edgeInd) = rootComplexGraph.Edges.Weight(edgeInd)./notConnectedEdgesValue;
        end
        bulkVector(newCellInd) = 1;
        rootComplexList = [rootComplexList newCellInd];
        newCellsInd = [newCellsInd newCellInd];
        
        
    end
    if(numel(newCellInd) == 0)
        t = max(amountNewCells, numel(outerborderInd));
        pressurePointsInd = outerborderInd(1:t);   
    else
        k = find(freeCellsInd == newCellInd);
        pressurePointsInd = freeCellsInd(1:k-1);   
    end
    
    %% Cells which couldnt grow are converted to pressure points
    
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
    %weiß nicht ob nötig:
    rootPressureDistributionVector(pressurePointsInd) = 1;
    %todo
    %rootPressureDistributionVector(rootPressureDistributionVector > 0) = normalize(rootPressureDistributionVector(rootPressureDistributionVector > 0),'range',[0 1]);
    
    
    optimalFutureRootComplexList = union(rootComplexList, pressurePointsInd);
    
    
    




%pressureCells werden zu potentiellen neuen Zellen
% da wo rootPressuredistribution == 1 ,da neue Zellen
end