function [rootComplexGraph, bulkVector, rootVector, mucilageVector, rootComplexParticleList, rootPressureDistributionVector] = ...
rootMucilageComplexGrowing(g, parameters, rootComplexGraph, bulkVector, rootVector, mucilageVector, rootComplexParticleList,  rootSourceCell, rootNr)

    %% Calculate Growing
    
%     currentAmountRootComplexCells = numel(rootComplexParticleList{rootNr});
%     
%     newAmountRootComplexCells = ceil(currentAmountRootCells * (parameters.rootGrowingRate + );
%     
    
    currentAmountRootCells = sum(rootVector, 'all');
    currentAmountMucilageCells = sum(mucilageVector, 'all');
    
    
    newAmountRootCells = ceil(currentAmountRootCells * parameters.rootGrowingRate);
    diffAmountRootCells = currentAmountRootCells - newAmountRootCells;
    
    newAmountMucilageCells = ceil(currentAmountMucilageCells * parameters.mucilageGrowingRate);
    diffAmountMucilageCells = currentAmountMucilageCells - newAmountMucilageCells;
    
    amountNewCells = diffAmountRootCells + diffAmountMucilageCells;
    
    %% Calculate possible Cells for Growing
    
    freeCellsInd = find((rootVector + mucilageVector) ~= 1);
    
    [TR,d] = shortestpathtree(rootComplexGraph,freeCellsInd,rootSourceCell);
    
    [sortedd, I] = sort(d);
    
    freeCellsInd = freeCellsInd(I);
    
    %
    i = 1;
    newCellsInd = [];
    while(amountNewCells >= numel(newCellsInd))
            st = stencil(g.NX,g.NX,freeCellsInd(i),1);
            %look if continuous growing without holes or skipped
            %connections
            isCellNeighboredRoot = sum(ismember(st(2:end),rootComplexParticleList{rootNr}),'all') > 0;
            if(isCellNeighboredRoot == 0)
                continue;
                %fprintf('geht nicht')
            end

            if(bulkVector(freeCellsInd(i)) == 0 && isCellNeighboredRoot)
                neighInd = neighbors(rootComplexGraph, freeCellsInd(i));    
                rootBorderCellsInd = intersect(neighInd,rootComplexParticleList{rootNr});

                h = freeCellsInd(i) * ones(size(rootBorderCellsInd,1),1);
                edgeInd = findedge(rootComplexGraph, rootBorderCellsInd, h);

                %nur van neumman weiterwachsen -> man dürfte nur Kanten mit 1 teilen
                %TODO sicher stellen dass überall richtige Werte
                if(rootComplexGraph.Edges.Weight(edgeInd) >= sqrt(2) * 1.1)
                    rootComplexGraph.Edges.Weight(edgeInd) = rootComplexGraph.Edges.Weight(edgeInd)./notConnectedEdgesValue;
                end
                bulkVector(freeCellsInd(i)) = 1;
                rootComplexParticleList{rootNr} = [rootComplexParticleList{rootNr} freeCellsInd(i)];
                newCellsInd = [newCellsInd freeCellsInd(i)];
            else
                %fprintf('Cell cant grow \n');
            end
            i = i+1;      
    end 
    %% Distinction between mucilage and root
    %totalAmountCells = numel(rootComplexParticleList{rootNr});
    
    
    
    
    
    
    %% Cells which couldnt grow are converted to pressure points
    visitedFreeCells = 1:i;
    growedFreeCells = (find(freeCellsInd == newCellsInd));
    pressurePointsInd = freeCellsInd(visitedFreeCells ~= growedFreeCells);
    
    
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
    pressureArea(rootComplexParticleList{rootNr}) = 0;
    rootPressureDistributionVector = rootPressureDistributionVector.*pressureArea;
    
    
    %todo
    %rootPressureDistributionVector(rootPressureDistributionVector > 0) = normalize(rootPressureDistributionVector(rootPressureDistributionVector > 0),'range',[0 1]);
    
    
    
    
    
    




%pressureCells werden zu potentiellen neuen Zellen
% da wo rootPressuredistribution == 1 ,da neue Zellen
end