function [rootComplexGraph, bulkVector, rootComplexList, rootPressureDistributionVector ,...
    bulkTypeVector, particleTypeVector, POMVector, POMconcVector, POMageVector, ...
    concAgent, concPOMAgent, POMagentAge,MucilageagentAge, edgeChargeVector, reactiveSurfaceVector, mucilageSurfaceVector,mucilageVector,...
    POMParticleList, solidParticleList, newCellsInd] = rootMucilageComplexGrowing...
    (g, rootComplexGraph, bulkVector, rootComplexVector,  rootComplexList, requieredAmountNewCells,...
    bulkTypeVector, particleTypeVector, POMVector, POMconcVector, POMageVector, ...
    concAgent, concPOMAgent, POMagentAge,MucilageagentAge, edgeChargeVector, reactiveSurfaceVector, mucilageSurfaceVector,mucilageVector,...
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
    
    freeSpaceStillPossible  = 1;
    index = 0;
    newCellGrowed = 1;
    cellOfInterestIndVector = [];
    while(freeSpaceStillPossible && numel(newCellsInd) < requieredAmountNewCells)
        if(newCellGrowed)
            %fprintf('newCelllGrowed \n')
            newCellGrowed = 0;
         %ahhh euklidisch
            freeCellsInd = find((rootComplexVector) ~= 1);

            [TR,d] = shortestpathtree(rootComplexGraph,freeCellsInd,rootSourceCell);

            [sortedd, I] = sort(d);
            freeCellsInd = freeCellsInd(I);

            outerborder = sortedd < notConnectedEdgesValue * 1.1;
            %TODO hier vielleicht nochmal eukliduscher Abstand berechnne
            outerborderInd = freeCellsInd(outerborder);%heeeeeeree 
            
            
            v = g.V0T(rootSourceCell,:);
            coordSource = mean(g.coordV(v,:),1);
            
            d = zeros(numel(outerborderInd),1);
            for i = 1:numel(outerborderInd)
                v = g.V0T(outerborderInd(i),:);
                coord = mean(g.coordV(v,:),1);
                X = [coord; coordSource];
                d(i) = pdist(X,'euclidean');
            end
            [sortedd, I] = sort(d);
            outerborderInd = outerborderInd(I);
            %d = pdist(X,'euclidean');
            freeSpots = find((bulkVector(outerborderInd) == 0));
            %cellOfInterestIndVector = [];%vielleicht auch hier
            if(numel(outerborderInd) == 0)
                break;
            end
        end
        index = index +1;
        %fprintf('index %d \n',index)
        if(index < numel(outerborderInd))
            freeSpaceStillPossible = 1;  
        else
            freeSpaceStillPossible = 0;
        end
        cellOfInterestInd = outerborderInd(index);
        cellOfInterestIndVector = [cellOfInterestIndVector, cellOfInterestInd];
        visitedNTimes = sum(cellOfInterestIndVector == cellOfInterestInd);
        if(visitedNTimes > 3)
            continue;
        end
        if(numel(newCellsInd) == ceil(requieredAmountNewCells/2)||numel(newCellsInd) == ceil(requieredAmountNewCells/4)||numel(newCellsInd) == ceil(requieredAmountNewCells/4 *3))
            cellOfInterestIndVector = [];
        end
        if(bulkVector(cellOfInterestInd) == 1)
              %find Connected Component
                bulkVector_bw = (reshape(bulkVector, [N N]));
                CC = bwconncomp(bulkVector_bw, 4);
                oo = cellfun(@(m) m == cellOfInterestInd,CC.PixelIdxList, 'UniformOutput', false);
                Match1 = cellfun(@sum, oo);
                connectedCompInd = find(Match1 == 1);
                
                %theoretisch müsste man weit wege Particle zu erst, letzter
                %der direkt an Wurzel
                %Find Solid ParticlesInd in this connected Area
                if(~isempty(connectedCompInd))
                    oo = cellfun(@(x) ismember(reshape(x',1,[]),CC.PixelIdxList{connectedCompInd}), solidParticleList, 'UniformOutput', false);
                    Match1 = cellfun(@sum, oo);
                    solidParticle = find(Match1 >= 1);
                end

                
                %Find POM ParticlesInd in this connected Area
                if(~isempty(connectedCompInd))
                    oo = cellfun(@(x) ismember(reshape(x',1,[]),CC.PixelIdxList{connectedCompInd}), POMParticleList, 'UniformOutput', false);
                    Match1 = cellfun(@sum, oo);
                    POMParticle = find(Match1 >= 1);
                end
                
                

                rootPressureDistributionVector = calculatePressureDistribution(g,cellOfInterestInd,bulkVector,rootComplexList);
                
                %rootPressureDistributionVector(:) = 0;
                rootPressureDistributionVector(freeSpots) = 1;
                if(numel(freeSpots) < requieredAmountNewCells)
					e = min(requieredAmountNewCells,numel(outerborderInd));
					rootPressureDistributionVector(outerborderInd(1:e)) = 1;
                end
                
                
                %solid Particle
                if(~isempty(solidParticle))
                    isTriedToPushAway = 0;%ismember(solidParticle, visitedSolidParticles);
                    if(~isTriedToPushAway)
                        for p = 1:numel(solidParticle)
                            particleSize = length( solidParticleList{ solidParticle(p) } ); 

                            fileID =0;
                            NZd = g.NX;
                            
                            bigParticleStencilLayers_individual = min(5, ceil(20/(particleSize)^0.5));
                            meanPressure = mean(rootPressureDistributionVector(solidParticleList{ solidParticle }));
                            bigParticleStencilLayers_individual = max(ceil(meanPressure * 5), bigParticleStencilLayers_individual);
                            
                            [bulkVector,bulkTypeVector, particleTypeVector, POMVector, POMconcVector, POMageVector,...
                                concAgent, concPOMAgent, POMagentAge,MucilageagentAge, edgeChargeVector, reactiveSurfaceVector,mucilageSurfaceVector,...
                                nextMucilageBorderVector, rootPressureDistributionVector,...
                                solidParticleList{ solidParticle(p) },~] = moveParticles( particleSize, bigParticleStencilLayers_individual,...
                                g, bulkVector, bulkTypeVector, particleTypeVector, POMVector, POMconcVector, POMageVector,...
                                concAgent, concPOMAgent, POMagentAge,MucilageagentAge, edgeChargeVector, reactiveSurfaceVector,mucilageSurfaceVector, ...
                                nextMucilageBorderVector, rootPressureDistributionVector, ...
                                NZd , fileID ,solidParticleList{ solidParticle(p) },0,4,0, attraction_type, 1);  

                           
                        end
                       erfolg = bulkVector(cellOfInterestInd) == 0;
                       if(erfolg == 1)
                           %fprintf('solid erfolg j %d \n',index)                      
                       else
                           %visitedSolidParticles = [visitedSolidParticles solidParticle(p)];
                       end
                        
                    end
                end                         
                %POM Particle
                if(~isempty(POMParticle))
                    isTriedToPushAway = 0;%ismember(POMParticle, visitedPOMParticles);
                    if(~isTriedToPushAway)
                        for p = 1:numel(POMParticle)
                            particleSize = length( POMParticleList{ POMParticle(p) } ); 

                            fileID =0;
                            NZd = g.NX;
                            
                            bigParticleStencilLayers_individual = min(5, ceil(20/(particleSize)^0.5));
                            meanPressure = mean(rootPressureDistributionVector(solidParticleList{ solidParticle(p) }));
                            bigParticleStencilLayers_individual = max(ceil(meanPressure * 5), bigParticleStencilLayers_individual);
                            
                            [bulkVector,bulkTypeVector, particleTypeVector, POMVector, POMconcVector, POMageVector,...
                                concAgent, concPOMAgent, POMagentAge,MucilageagentAge, edgeChargeVector, reactiveSurfaceVector,mucilageSurfaceVector,...
                                nextMucilageBorderVector, rootPressureDistributionVector,...
                                POMParticleList{ POMParticle(p) },~] = moveParticles( particleSize, bigParticleStencilLayers_individual,...
                                g, bulkVector, bulkTypeVector, particleTypeVector, POMVector, POMconcVector, POMageVector,...
                                concAgent, concPOMAgent, POMagentAge,MucilageagentAge, edgeChargeVector, reactiveSurfaceVector,mucilageSurfaceVector, ...
                                nextMucilageBorderVector, rootPressureDistributionVector, ...
                                NZd , fileID ,POMParticleList{ POMParticle(p) },0,4,0, attraction_type, 1);  
                        end
                       erfolg = bulkVector( cellOfInterestInd) == 0;
                       if(erfolg == 1)
                          % fprintf('Pom erfolg j %d \n',index)
                        else
                          % visitedPOMParticles = [visitedPOMParticles POMParticle(p)];
                       end
                    end   
                end %POM Particle
                
        end        
               %index = 0;%oder vielleicht index = index -1
        if(bulkVector(cellOfInterestInd) == 0)
             %Markieren
            fprintf('grow %d \n',index)
            newCellInd = cellOfInterestInd;
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
            rootComplexVector(newCellInd) = 1;
            rootComplexList = [rootComplexList newCellInd];
            newCellsInd = [newCellsInd newCellInd];
            
            
            index = 0;
            newCellGrowed = 1;
        else
            %klappt nicht 
        end 
    end
   pressurePointsInd = setdiff(cellOfInterestIndVector,newCellsInd);
%    if(numel(newCellsInd) < requieredAmountNewCells)
%        diff = requieredAmountNewCells - numel(newCellsInd);
%        pressurePointsInd = outerborderInd(:);
%    else
%        pressurePointsInd = [];
%    end

      
    %% Cells which couldnt grow are converted to pressure points
    rootPressureDistributionVector = calculatePressureDistribution(g,pressurePointsInd,bulkVector,rootComplexList);
end


