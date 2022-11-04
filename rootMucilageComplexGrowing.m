function [rootComplexGraph, bulkVector, rootComplexList, rootPressureDistributionVector , optimalFutureRootComplexList,...
    bulkTypeVector, particleTypeVector, POMVector, POMconcVector, POMageVector, ...
    concAgent, concPOMAgent, POMagentAge, edgeChargeVector, reactiveSurfaceVector, mucilageVector,...
    POMParticleList, solidParticleList] = rootMucilageComplexGrowing...
    (g, rootComplexGraph, bulkVector, rootComplexVector,  rootComplexList, amountNewCells,...
    bulkTypeVector, particleTypeVector, POMVector, POMconcVector, POMageVector, ...
    concAgent, concPOMAgent, POMagentAge, edgeChargeVector, reactiveSurfaceVector, mucilageVector,...
    POMParticleList, solidParticleList,...
     attraction_type)

    N = g.NX;
    rootSourceCell = rootComplexList(1);
    notConnectedEdgesValue = N* N * 2;
    newCellsInd = [];
    %das ielleicht noch besser machen wie im haupteil
    nextMucilageBorderVector = mucilageVector;
    rootPressureDistributionVector = 0 * ones(g.numT, 1);
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
        
        k = find((bulkVector(outerborderInd) == 0),1,'first');
        pushaway = 0;
        if(~isempty(k))
            if(k > 1)
            pushaway = 1;
            end
        else
            pushaway = 1;
            k = numel(outerborderInd);
        end
        
        
        if(pushaway)
           
            occupiedCellsInd = outerborderInd(bulkVector(outerborderInd) == 1);
            for j = 1:k-1
                
                rootPressureDistributionVector(occupiedCellsInd(j)) = 1;
                oo = cellfun(@(x) x==occupiedCellsInd(j), solidParticleList, 'UniformOutput', 0);
                Match1 = cellfun(@sum, oo);
                solidParticle = find(Match1 == 1);

                if(~isempty(solidParticle))
                    particleSize = length( solidParticleList{ solidParticle } ); 

                    fileID =0;
                    NZd = g.NX;
                    bigParticleStencilLayers_individual = min(5, ceil(20/(particleSize)^0.5));
                    [bulkVector,bulkTypeVector, particleTypeVector, POMVector, POMconcVector, POMageVector,...
                        concAgent, concPOMAgent, POMagentAge, edgeChargeVector, reactiveSurfaceVector,...
                        nextMucilageBorderVector, rootPressureDistributionVector,...
                        solidParticleList{ solidParticle },~] = moveParticles( particleSize, bigParticleStencilLayers_individual,...
                        g, bulkVector, bulkTypeVector, particleTypeVector, POMVector, POMconcVector, POMageVector,...
                        concAgent, concPOMAgent, POMagentAge, edgeChargeVector, reactiveSurfaceVector, ...
                        nextMucilageBorderVector, rootPressureDistributionVector, ...
                        NZd , fileID ,solidParticleList{ solidParticle },0,4,0, attraction_type);  

                   erfolg = bulkVector( occupiedCellsInd(j)) == 0;
                   if(erfolg == 1)
                       break;
                   end
                end
            end
        end
        
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
    rootPressureDistributionVector(:) = 0;
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