function [rootComplexGraph, bulkVector, rootComplexList, rootPressureDistributionVector ,...
    bulkTypeVector, particleTypeVector, POMVector, POMconcVector, POMageVector, ...
    concAgent, concPOMAgent, POMagentAge,MucilageagentAge, edgeChargeVector, reactiveSurfaceVector, mucilageSurfaceVector,mucilageVector,...
    POMParticleList, solidParticleList, newCellsInd] = rootMucilageComplexGrowing_3...
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
    
    maximum_width = 1;%35;

    %% Calculate possible Cells for Growing
    
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
                if(abs(coord(2) - coordSource(2)) > maximum_width)
                    %fprintf("Weg")
                    d(i) = inf;
                end
            end
            
            [sortedd, I] = sort(d);
            %TODO sorted OuterborderI
            outerborderInd = outerborderInd(I);
            %d = pdist(X,'euclidean');
            freeSpots = find((bulkVector(outerborderInd) == 0));
            %cellOfInterestIndVector = [];%vielleicht auch hier
            if(numel(outerborderInd) == 0  || sum(d == inf) == numel(outerborderInd) )%
                break;
            end
        end
        index = index +1;
        %fprintf('index %d \n',index)
        %if(index < numel(outerborderInd))
            freeSpaceStillPossible = 1;  
        %else
            %freeSpaceStillPossible = 0;
        %end
         if(all(d == inf))
            break;
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
        
        
        
        
        %----MOVE--Start------------
        

        
        %MOoooooooooove         
        if(bulkVector(cellOfInterestInd) == 1)

            [particleList, particleContent] = particleInfoTUM(bulkVector-(rootComplexVector), solidParticleList, POMParticleList);
             oo = cellfun(@(m) m == cellOfInterestInd,particleList, 'UniformOutput', false);
             Match1 = cellfun(@sum, oo);
             connectedCompInd = find(Match1 == 1);

               aggInd = particleList{connectedCompInd};
              if(~isempty(connectedCompInd))
                        oo = cellfun(@(x) ismember(reshape(x',1,[]),particleList{connectedCompInd}), solidParticleList, 'UniformOutput', false);
                        Match1 = cellfun(@sum, oo);
                        connSolidParticle = find(Match1 >= 1);
                        
                        oo = cellfun(@(m) m == cellOfInterestInd,solidParticleList, 'UniformOutput', false);
                        Match1 = cellfun(@sum, oo);
                        refIndPart = find(Match1 == 1);
                        
              end
              
                %Find POM ParticlesInd in this connected Area
                if(~isempty(connectedCompInd))
                    oo = cellfun(@(x) ismember(reshape(x',1,[]),particleList{connectedCompInd}), POMParticleList, 'UniformOutput', false);
                    Match1 = cellfun(@sum, oo);
                    connPOMParticle = find(Match1 >= 1);
                    
                        oo = cellfun(@(m) m == cellOfInterestInd,POMParticleList, 'UniformOutput', false);
                        Match1 = cellfun(@sum, oo);
                        refIndPOM = find(Match1 == 1);
                end 
                
                
              bulkVector_agg(:) = bulkVector;
              bulkVector_agg(:) = 0;
              bulkVector_agg(aggInd) = 1;
              dMap = distanceMap(bulkVector_agg, cellOfInterestInd);
              %dMap(dMap ~=0) = inf;
              dMap_geo = reshape(dMap, [N,N]);  
              %imshow(dMap_geo, [0,100])
              
              
              agg = particleList{connectedCompInd};
              
               minDsolid = zeros(numel(connSolidParticle),1);
              for i = 1:numel(connSolidParticle)
                  solid = solidParticleList{connSolidParticle(i)};
                  minDsolid(i) = min(dMap(solid));
              end
              minDPOM = zeros(numel(connPOMParticle),1);
              for i = 1:numel(connPOMParticle)
                  pom = POMParticleList{connPOMParticle(i)};
                  minDPOM(i) = min(dMap(pom));
              end
              sortedSizes = unique(sort([minDPOM;minDsolid]));
              
              for i = 1:numel(sortedSizes)
                  bulkVector_agg = bulkVector;
                  bulkVector_agg(:) = 0;
                  bulkVector_agg(aggInd) = 1;
                  try
                  particles = connSolidParticle(find(minDsolid >= sortedSizes(i)));
                  catch 
                      find(minDsolid >= sortedSizes(i))
                  end
                  poms = connPOMParticle(find(minDPOM >= sortedSizes(i)));
                  if(~isempty(particles))
                    for j = 1:numel(particles)
                    bulkVector_agg(solidParticleList{particles(j)}) = 0;
                    end
                  end
                  
                  if(~isempty(poms))
                    for j = 1:numel(poms)
                    bulkVector_agg(POMParticleList{poms(j)}) = 0;
                    end
                  end
                  
                  subb_agg(i) = {bulkVector_agg};
                  
                  %bulkVector_agg_geo = reshape(bulkVector_agg, [N,N]);  
                  %imshow(bulkVector_agg_geo)
                  
                  
              end
              
             
              
             
                rootPressureDistributionVector = calculatePressureDistribution(g,cellOfInterestInd,bulkVector,rootComplexList);
                
                %rootPressureDistributionVector(:) = 0;
                rootPressureDistributionVector(freeSpots) = 1;
                if(numel(freeSpots) < requieredAmountNewCells)
					e = min(requieredAmountNewCells,numel(outerborderInd));
					rootPressureDistributionVector(outerborderInd(1:e)) = 1;
                end 
              

				[particleList, particleContent] = particleInfoTUM(subb_agg{i}, solidParticleList, POMParticleList);
				for particle = 1 : length( particleList )
					particleSize = length( particleList{ particle } ); 
					fileID =0;
					NZd = g.NX;
								
					if(size(particleContent{particle},1)<2)% kein Verbund
						continue
					end
					if(size(particleContent{particle},1) == 10)% contains root
						continue
					end
					test_ind = particleList{particle}(1);

					bigParticleStencilLayers_individual = min(5, ceil(20/(particleSize)^0.5));
					meanPressure = mean(rootPressureDistributionVector( particleList{ particle }));
					bigParticleStencilLayers_individual = max(ceil(meanPressure * 5), bigParticleStencilLayers_individual);

					bigParticleStencilLayers_individual = 5;
					if particleSize > 20000
						bigParticleStencilLayers_individual = 0;
					end
					if bigParticleStencilLayers_individual == 0
						continue
					end
					%bulkMucilageVector = bulkVector + mucilageVector
								
					 bulkVector_before = bulkVector;    
					[bulkVector,bulkTypeVector, particleTypeVector, POMVector, POMconcVector, POMageVector, ...
						concAgent, concPOMAgent, POMagentAge,MucilageagentAge, edgeChargeVector, reactiveSurfaceVector,mucilageSurfaceVector,...
						mucilageVector, rootPressureDistributionVector,...
						changedList,~] = ...
						moveParticles( particleSize, bigParticleStencilLayers_individual, g, bulkVector, bulkTypeVector, ... 
						particleTypeVector, POMVector, POMconcVector, POMageVector, concAgent, ...
						concPOMAgent, POMagentAge, MucilageagentAge, edgeChargeVector, reactiveSurfaceVector, mucilageSurfaceVector, ...
						mucilageVector, rootPressureDistributionVector, ...    
						NZd , fileID ,particleList{ particle },0,4,0, attraction_type, 1); 
				   %anser = sum(bulkVector_before ~= bulkVector)
				   % bulkVector = bulkMucilageVector - mucilageVector;

					if test_ind ~= changedList(1) % Listen müssen angepasst werden
						sten_ind = find(stencil(g.NX,g.NX,test_ind,bigParticleStencilLayers_individual)==changedList(1));
						for part = 1: size(particleContent{particle},1)  
							switch particleContent{particle}(part,1)
								case 1 
									sten = stencil(g.NX,g.NX,solidParticleList{ particleContent{particle}(part,2)},bigParticleStencilLayers_individual);
									solidParticleList{particleContent{particle}(part,2)}=sten(:,sten_ind)';
								case 2 
									sten = stencil(g.NX,g.NX,POMParticleList{ particleContent{particle}(part,2)},bigParticleStencilLayers_individual);
									POMParticleList{particleContent{particle}(part,2)}=sten(:,sten_ind)';
							end
						end
						break;
				   
					end 
					break;

				end	


                solidParticle = connSolidParticle;
                %solid Particle
                if(~isempty(solidParticle))
                    isTriedToPushAway = 0;%ismember(solidParticle, visitedSolidParticles);
                    if(~isTriedToPushAway)
                        for p = 1:numel(solidParticle)
                            particleSize = length( solidParticleList{ solidParticle(p) } ); 

                            fileID =0;
                            NZd = g.NX;
                            
                            bigParticleStencilLayers_individual = min(5, ceil(20/(particleSize)^0.5));
                            meanPressure = mean(rootPressureDistributionVector(solidParticleList{ solidParticle(p) }));
                            bigParticleStencilLayers_individual = max(ceil(meanPressure * 5), bigParticleStencilLayers_individual);
                            bigParticleStencilLayers_individual  = 10;
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
							%break;						   
                       else
                           %visitedSolidParticles = [visitedSolidParticles solidParticle(p)];
                       end
                        
                    end
                end                         
                %POM Particle
                POMParticle = connPOMParticle;
                if(~isempty(POMParticle))
                    isTriedToPushAway = 0;%ismember(POMParticle, visitedPOMParticles);
                    if(~isTriedToPushAway)
                        for p = 1:numel(POMParticle)
                            particleSize = length( POMParticleList{ POMParticle(p) } ); 

                            fileID =0;
                            NZd = g.NX;
                            
                            bigParticleStencilLayers_individual = min(5, ceil(20/(particleSize)^0.5));
                            meanPressure = mean(rootPressureDistributionVector(POMParticleList{ POMParticle(p) }));
                            bigParticleStencilLayers_individual = max(ceil(meanPressure * 5), bigParticleStencilLayers_individual);
                            bigParticleStencilLayers_individual  = 10;
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
						   %break;
                        else
                          % visitedPOMParticles = [visitedPOMParticles POMParticle(p)];
						 
                       end
                    end   
                end %POM Particle


                
                
        end 
        
        
        
        
        
        %----MOVE---END--------------
        
        
        if(bulkVector(cellOfInterestInd) == 0)
             %Markieren
            fprintf('grow %d number %d \n',index, numel(newCellsInd))
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
fprintf("Ende Grow \n")
end


