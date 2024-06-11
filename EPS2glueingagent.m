 function [bulkVector, EPSconcVector, concPOMAgent, edgeChargeVector] = EPS2glueingagent(g, parameters, bulkVector, EPSVector, EPSconcVector,concPOMAgent, edgeChargeVector)
    
    EPSInd = find(EPSVector == 1);
    EPSParticleList = cell(numel(EPSInd),1);
    for i = 1:numel(EPSInd)
        EPSParticleList{i} = EPSInd(i);
    end
    
    
    % here concPOMAgent is considered as an edgeChargeVector
    EPSsolidEdgeList = calculatePOMsolidEdgeList(g, bulkVector, EPSVector, EPSParticleList);
   
    
   for i = 1 : length(EPSParticleList)
        if ~isempty(EPSsolidEdgeList{i})
        concPOMAgent(EPSsolidEdgeList{i}) = concPOMAgent(EPSsolidEdgeList{i}) + 0.025;    
        end
    end
    
    
    % spreading
    for i = 1 : length(EPSsolidEdgeList)
        for edge = 1 : length(EPSsolidEdgeList{i})
           if concPOMAgent(EPSsolidEdgeList{i}(edge)) > parameters.POMagentMax
               excessPOM = concPOMAgent(EPSsolidEdgeList{i}(edge)) - parameters.POMagentMax;
               concPOMAgent(EPSsolidEdgeList{i}(edge)) = parameters.POMagentMax;
               visitedVertices = g.V0CE(EPSsolidEdgeList{i}(edge),:);
               visitedEdges = EPSsolidEdgeList{i}(edge);
               [edgeHelper, ~] = find(ismember(g.V0CE, visitedVertices));
               edgeHelper = edgeHelper(~ismember(edgeHelper,visitedEdges));
               edgeCandidates = [];
               for indEdge = 1 : length(edgeHelper)
                  solidInd = mod(edgeHelper(indEdge), g.NX * g.NX);
                  if ((bulkVector(solidInd) == 1) &&  (POMVector(solidInd) == 0))
                      possibleNeighbors = stencil( g.NX, g.NX, solidInd, 1); 
                      possibleNeighbors = possibleNeighbors(2:end); 
                      [~,edgeDirection] = find(g.CE0T(solidInd,:)==edgeHelper(indEdge));
                      for neigh = 1 : 4
                        if ( (neigh == 1) && (edgeDirection == 1) )
                            if ((POMVector(possibleNeighbors(neigh)) == 1) || (bulkVector(possibleNeighbors(neigh)) == 0))
                                edgeCandidates = [edgeCandidates edgeHelper(indEdge)];
                            end
                        elseif ( (neigh == 2) && (edgeDirection == 4) )
                            if ((POMVector(possibleNeighbors(neigh)) == 1) || (bulkVector(possibleNeighbors(neigh)) == 0))
                                edgeCandidates = [edgeCandidates edgeHelper(indEdge)];
                            end
                        elseif ( (neigh == 3) && (edgeDirection == 3) )
                            if ((POMVector(possibleNeighbors(neigh)) == 1) || (bulkVector(possibleNeighbors(neigh)) == 0))
                                edgeCandidates = [edgeCandidates edgeHelper(indEdge)];
                            end
                        elseif ( (neigh == 4) && (edgeDirection == 2) )
                            if ((POMVector(possibleNeighbors(neigh)) == 1) || (bulkVector(possibleNeighbors(neigh)) == 0))
                                edgeCandidates = [edgeCandidates edgeHelper(indEdge)];
                            end
                        end
                      end
                  end
               end
             
               while ( (excessPOM > 0) &&  ~isempty(edgeCandidates) )
                   % get capacity
                   aim = edgeCandidates(1);
                   edgeCandidates = edgeCandidates(2:end);
                   visitedEdges = [visitedEdges aim];
                   
                   aimCapacity = max(parameters.POMagentMax - concPOMAgent(aim), 0);
                   
                   % edgeCandidate has enough capacity to take up the
                   % excessPOM
                   if(aimCapacity >= excessPOM)
                       concPOMAgent(aim) = concPOMAgent(aim) + excessPOM;
                       excessPOM = 0;
                   % aim can not take up all the excess POM    
                   else
                       if aimCapacity > 0
                          concPOMAgent(aim) = parameters.POMagentMax;
                       end
                       excessPOM = excessPOM - aimCapacity;
                       
                       verticesCandidates = g.V0CE(aim, :);
                       verticesCandidates = verticesCandidates(~ismember(verticesCandidates,visitedVertices));
                       % verticesCandidates should never be more than one,
                       % because every vertex comes from one that has
                       % already been visited
                       if ~isempty(verticesCandidates)
                           [edgeHelper, ~] = find(ismember(g.V0CE, verticesCandidates));
                           edgeHelper = edgeHelper(~ismember(edgeHelper,visitedEdges));
                           
                           for indEdge = 1 : length(edgeHelper)
                              solidInd = mod(edgeHelper(indEdge), g.NX * g.NX);
                              if ((bulkVector(solidInd) == 1) &&  (POMVector(solidInd) == 0))
                                  possibleNeighbors = stencil( g.NX, g.NX, solidInd, 1); 
                                  possibleNeighbors = possibleNeighbors(2:end); 
                                  [~,edgeDirection] = find(g.CE0T(solidInd,:)==edgeHelper(indEdge));
                                  for neigh = 1 : 4
                                    if ( (neigh == 1) && (edgeDirection == 1) )
                                        if ((POMVector(possibleNeighbors(neigh)) == 1) || (bulkVector(possibleNeighbors(neigh)) == 0))
                                            edgeCandidates = [edgeCandidates edgeHelper(indEdge)];
                                        end
                                    elseif ( (neigh == 2) && (edgeDirection == 4) )
                                        if ((POMVector(possibleNeighbors(neigh)) == 1) || (bulkVector(possibleNeighbors(neigh)) == 0))
                                            edgeCandidates = [edgeCandidates edgeHelper(indEdge)];
                                        end
                                    elseif ( (neigh == 3) && (edgeDirection == 3) )
                                        if ((POMVector(possibleNeighbors(neigh)) == 1) || (bulkVector(possibleNeighbors(neigh)) == 0))
                                            edgeCandidates = [edgeCandidates edgeHelper(indEdge)];
                                        end
                                    elseif ( (neigh == 4) && (edgeDirection == 2) )
                                        if ((POMVector(possibleNeighbors(neigh)) == 1) || (bulkVector(possibleNeighbors(neigh)) == 0))
                                            edgeCandidates = [edgeCandidates edgeHelper(indEdge)];
                                        end
                                    end
                                  end
                              end
                           end
                       end
                   end
                                   
                   
               end
               if (excessPOM > 0)
                    fprintf('POM agent could not be spread!')
                    sumExcessPOM = sumExcessPOM + excessPOM;
                    excessPOM = 0;
               end
%                assert(excessPOM == 0, 'POM agent could not be spread!')
               
           end
        end
    end
    

    
    edgeChargeVector( concPOMAgent > parameters.POMagentMin ) = 1;



    
 
 end
