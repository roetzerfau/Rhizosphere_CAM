function [mucilageConcVector, mucilageVector, mucilageGraph ] = updateMucilage(g, p, bulkVector, outerRootBorderInd, mucilageConcVector, mucilageVector, mucilageGraph)
    
    %----------------------------------------------------------------
    %bei wegschieben, bekommt nächstgelende freie(kein Bulk, beliebige Konzetration) Zelle Con
    %freie Zelle = kein Bulk
    %hier nicht unbedingt keine Border celle
    mucilageConcVector_before = mucilageConcVector;
    mucilageAmount_before = sum(mucilageConcVector_before);
    
    %mucilageConcVector(bulkVector == 1) = 0;
    distributedMucilageInd = find(bulkVector ==1 & mucilageVector ==1);
    for i = 1:numel(distributedMucilageInd)
        mucilageConcVector(distributedMucilageInd(i))= 0;
        overshootValue = mucilageConcVector_before(distributedMucilageInd(i));
        nextFreeCellsInd = findNextFreeCells(bulkVector, mucilageConcVector, mucilageGraph,distributedMucilageInd(i), 5, 0);
        if(numel(nextFreeCellsInd) > 0)
            mucilageConcVector(nextFreeCellsInd(1)) = mucilageConcVector(nextFreeCellsInd(1)) + overshootValue;
        else
            %mucilageConcVector = mucilageConcVector_before;
        end
    end
    mucilageAmount_after = sum(mucilageConcVector);
    if(abs(mucilageAmount_before - mucilageAmount_after) > 0.0001)
          error('Falsch wegdrücken')
    end
    [mucilageVector, mucilageConcVector, mucilageGraph] = updateMucilageVector(p.minConcMucilage , mucilageVector, mucilageConcVector, mucilageGraph);
    %------------------------------------------------------------------------------------------
    %%Wachsen von Zellen am Rand
    mucilageConcVector_before = mucilageConcVector;
    mucilageAmount_before = sum(mucilageConcVector_before);
    if(p.mucilageGrowing == 1)
        for i = 1:numel(outerRootBorderInd)
            if( bulkVector(outerRootBorderInd(i)) == 0 && mucilageConcVector(outerRootBorderInd(i)) <= 5 )
               mucilageConcVector(outerRootBorderInd(i)) = mucilageConcVector(outerRootBorderInd(i)) + 1.0;
            end
        end
        [mucilageVector, mucilageConcVector, mucilageGraph] = updateMucilageVector(p.minConcMucilage , mucilageVector, mucilageConcVector, mucilageGraph);
    end 
    mucilageAmount_after = sum(mucilageConcVector);
    if(abs(mucilageAmount_before - mucilageAmount_after) > numel(outerRootBorderInd))
          error('Falsch wegdrücken')
    end
    %------------------------------------------------------------------------------------------
    
    %(Cellen an Wurzel && freie Zelle) mmüssen Conc 1 haben
    %für jede Zelle wird nächstgelenden verbundene BorderZelle gesucht. 
    %von den so viel abzapfen bis wieder 1. 
%     mucilageConcVector_before = mucilageConcVector;
%     mucilageAmount_before = sum(mucilageConcVector_before);
%     
%     ind_borderfalsch = find((mucilageConcVector(outerRootBorderInd) < 1 + ~(bulkVektor >= 1)) == 2);
%     for i = 1: numel(ind_borderfalsch)
%         %hier nächste Border Zelle finden, egal welche Konzentration
%     end
%     
%     
%     mucilageAmount_after = sum(mucilageConcVector);
%     if(mucilageAmount_before ~= mucilageAmount_after)
%         error('Falsch Border')
%     end
%     
    
    %------------------------------------------------------------------------------------------
    
    %suche nach von Überschuss -> zelle nächste verbundene Zelle, aber kleiner als 1;
    mucilageConcVector_before = mucilageConcVector;
    mucilageAmount_before = sum(mucilageConcVector_before);
    overshootConcInd = find(mucilageConcVector > 1);
    for i = 1:numel(overshootConcInd)
        overshootValue = mucilageConcVector(overshootConcInd(i))-1;
        %next free cell muss verbunden sein, also border cell
        nextFreeCellsInd = findNextFreeCells(bulkVector, mucilageConcVector, mucilageGraph, overshootConcInd(i), 1, 1);
        index  = 0;
        while(overshootValue >=0)
            index = index + 1;
            if(numel(nextFreeCellsInd) < index)
                break;
            end
            conNFC = mucilageConcVector(nextFreeCellsInd(index));
            mucilageConcVector(nextFreeCellsInd(index)) = min (1, conNFC + overshootValue);
            diffConc = mucilageConcVector(nextFreeCellsInd(index)) - conNFC;
            overshootValue = overshootValue - diffConc;
            mucilageConcVector(overshootConcInd(i)) = 1 + overshootValue;
        end
    end
     mucilageAmount_after = sum(mucilageConcVector);
    if(abs(mucilageAmount_before - mucilageAmount_after) > 0.0001)
          error('Falsch wegdrücken')
    end
    [mucilageVector, mucilageConcVector, mucilageGraph] = updateMucilageVector(p.minConcMucilage , mucilageVector, mucilageConcVector, mucilageGraph);

     [mucilageConcVector] = calculateMucilageDecay(g, p.mucilageDecayRate, bulkVector, mucilageVector, mucilageConcVector);
     [mucilageVector, mucilageConcVector, mucilageGraph] = updateMucilageVector(p.minConcMucilage , mucilageVector, mucilageConcVector, mucilageGraph);
    
    %----------------------------------------------------------------
end

function ind = findNextFreeCells(bulkVector, mucilageConcVector, mucilageGraph, startCell, targetConc, isDirectConnected)
            global notConnectedEdgesValue;
            freeCellsInd = find((bulkVector + ~(mucilageConcVector < targetConc)) == 0);

            [TR,d] = shortestpathtree(mucilageGraph,freeCellsInd,startCell);

            [sortedd, I] = sort(d);
            freeCellsInd = freeCellsInd(I);
            %ind = freeCellsInd;
            
            %mzaa verbundn sein 
            %stilll connected, also entweder ist mucilage mi weniger
            %konzentration oder direkter nachbar
            if(isDirectConnected ==1)
                directConnected = sortedd < notConnectedEdgesValue * 1.1;  
                %TODO hier vielleicht nochmal eukliduscher Abstand berechnne
                directConnectedInd = freeCellsInd(directConnected);%heeeeeeree
                ind = directConnectedInd;
            else
                ind = freeCellsInd;
            end
            
end
function [mucilageVector, mucilageConcVector, mucilageGraph] = updateMucilageVector(minConcMucilage , mucilageVector, mucilageConcVector, mucilageGraph)
 
    % Remove
    indMucilagechanged = find( ( mucilageVector == 1 ) &  ( mucilageConcVector < minConcMucilage ) );
    if(~isempty(indMucilagechanged))
        for i = 1:numel(indMucilagechanged)
            oldCellInd = indMucilagechanged(i);
            [mucilageGraph, mucilageVector] = removeMucilageCell(mucilageGraph, mucilageVector,oldCellInd);
            mucilageConcVector(oldCellInd) = 0;
        end
    end
    % Add
    indMucilagechanged = find( ( mucilageVector == 0 ) &  ( mucilageConcVector > minConcMucilage));
    if(~isempty(indMucilagechanged))
        for i = 1:numel(indMucilagechanged)
            newCellInd = indMucilagechanged(i);
            [mucilageGraph, mucilageVector] = addMucilageCell(mucilageGraph, mucilageVector, newCellInd);
        end
    end
end
function [mucilageGraph, mucilageVector] = addMucilageCell(mucilageGraph, mucilageVector, newCellInd)
   
    global notConnectedEdgesValue;
    mucilageList = find(mucilageVector);
    neighInd = neighbors(mucilageGraph, newCellInd);    
    mucilageNeighCellsInd = intersect(neighInd,mucilageList);

    h = newCellInd * ones(size(mucilageNeighCellsInd,1),1);
    edgeInd = findedge(mucilageGraph, mucilageNeighCellsInd, h);

    if(mucilageGraph.Edges.Weight(edgeInd) >= sqrt(2) * 1.1)
        mucilageGraph.Edges.Weight(edgeInd) = mucilageGraph.Edges.Weight(edgeInd)./notConnectedEdgesValue;
    end
    mucilageVector(newCellInd) = 1;
    
end
function [mucilageGraph, mucilageVector] = removeMucilageCell(mucilageGraph, mucilageVector, oldCellInd)
    global notConnectedEdgesValue;


    mucilageList = find(mucilageVector);
    neighInd = neighbors(mucilageGraph, oldCellInd);    
    mucilageNeighCellsInd = intersect(neighInd,mucilageList);

    h = oldCellInd * ones(size(mucilageNeighCellsInd,1),1);
    edgeInd = findedge(mucilageGraph, mucilageNeighCellsInd, h);
    
    if(mucilageGraph.Edges.Weight(edgeInd) < sqrt(2) * 1.1)
        mucilageGraph.Edges.Weight(edgeInd) = mucilageGraph.Edges.Weight(edgeInd).* notConnectedEdgesValue;
    end
    mucilageVector(oldCellInd) = 0;
end