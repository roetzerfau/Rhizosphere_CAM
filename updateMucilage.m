function [] = updateMucilage(bulkVector, outerRootBorderInd, mucilageConcVector, mucilageGraph)
    %jede Zelle wächst 
    % grwoing rate * 1 wenn an root
    % find outer Border
    % wenn conc nicht 1, dann wachsen
    % find nächstgelegende zellenstelle/ wenn weg größer als wert nicht
    % verbunden dann kein wachstum 
    % 
    %insgesamte Mucilage Conc
    
    % überschüssige konzetration in nächstgelende freie Zelle rein
    %wurzel rückgang: so viel neu reinpumpen
    %als mucilage markierte Feld müssen immer conc haben
    
    %von größten überschuss abwärts 
    
    %cellen an wurzel/mucilage sind mucilage_potenzial gekennzeichnet 
    %geringste als mucliage markierte Zelle bekommt entweder überschuss oder unterschuss am nächstgelende
    %celle 

    
    
    %----------------------------------------------------------------
    %bei wegschieben, bekommt nächstgelende freie(kein Bulk, beliebige Konzetration) Zelle Con
    %freie Zelle = kein Bulk
    %hier nicht unbedingt keine Border celle
    mucilageConcVector_before = mucilageConcVector;
    mucilageAmount_before = sum(mucilageConcVector_before);
    
    mucilageConcVector(bulkVector == 1) = 0;
    distributedMucilageInd = find(mucilageConcVector_before ~=mucilageConcVector);
    for i = 1:numel(distributedMucilageInd)
        overshootValue = mucilageConcVector_before(distributedMucilageInd(i));
        nextFreeCellsInd = findNextFreeCells(bulkVector, mucilageConcVector, mucilageGraph);
        if(numel(FreeCellsInd) > 0)
            mucilageConcVector(nextFreeCellsInd(1)) = mucilageConcVector(nextFreeCellsInd(1)) + overshootValue;
        end
    end
    mucilageAmount_after = sum(mucilageConcVector);

    if(mucilageAmount_before ~= mucilageAmount_after)
          error('Falsch wegdrücken')
    end
    mucilageVector = updateMucilageVector();
    
    %%Wachsen von Zellen am Rand
    
    
    
    %(Cellen an Wurzel && freie Zelle) mmüssen Conc 1 haben
    %für jede Zelle wird nächstgelenden verbundene BorderZelle gesucht. 
    %von den so viel abzapfen bis wieder 1. 
    mucilageConcVector_before = mucilageConcVector;
    mucilageAmount_before = sum(mucilageConcVector_before);
    
    ind_borderfalsch = find((mucilageConcVector(outerRootBorderInd) < 1 + ~(bulkVektor >= 1)) == 2);
    for i = 1: numel(ind_borderfalsch)
        %hier nächste Border Zelle finden, egal welche Konzentration
    end
    
    
    mucilageAmount_after = sum(mucilageConcVector);
    if(mucilageAmount_before ~= mucilageAmount_after)
        error('Falsch Border')
    end
    
    
   
    
    %suche nach von Überschuss -> zelle nächste verbundene Zelle, aber kleiner als 1;
    overshootConcInd = find(mucilageConcVector > 1);
    for i = 1:numel(overshootConcInd)
        overshootValue = mucilageConcVector(overshootConcInd(i))-1;
        %next free cell muss verbunden sein, also border cell
        nextFreeCellsInd = findNextFreeCells(bulkVector, mucilageConcVector, mucilageGraph, overshootConcInd(i));
        index  = 0;
        while(overshootValue ~=0)
            index = index + 1;
            conNFC = mucilageConcVector(nextFreeCellsInd(index));
            mucilageConcVector(nextFreeCellsInd(index)) = conNFC + min (1, conNFC + overshootValue);
            diffConc = conNFC - mucilageConcVector(nextFreeCellsInd(index));
            overshootValue = overshootValue - diffConc;
        end
    end
    
    %man kann durch conc 0.001 was markieren
    
    
    

    %----------------------------------------------------------------


    
    
    updateCellsAmount = sum(mucilageConcVector_before - mucilageConcVector);
    
    
    %outerborder
    freeCellsInd = find((rootComplexVector) ~= 1);
    [TR,d] = shortestpathtree(rootComplexGraph,freeCellsInd,rootSourceCell);
    [sortedd, I] = sort(d);
    freeCellsInd = freeCellsInd(I);
    
    if(requiredAmountChangeCells < 0)
        
    elseif(requiredAmountChangeCells > 0)
    else
    end
%     outerborder = sortedd < notConnectedEdgesValue * 1.1;
%     outerborderInd = freeCellsInd(outerborder);%heeeeeeree 
    
    
    
    
    %das zählen und auf required Amount of new clless draufrechnen
    
    %wieder kompakt machen
    
    % und dann mit bilanz weiterarbeiten
    %im Endeffekt einfach nur plus oder minus cellen
    
    %kenzentrationen verschieben
end

function ind = findNextFreeCells(bulkVector, mucilageConcVector, mucilageGraph, startCell, targetConc)
            notConnectedEdgesValue = N* N * 2;
            freeCellsInd = find((bulkVector + (mucilageConcVector >= targetConc)) ~= 1);

            [TR,d] = shortestpathtree(mucilageGraph,freeCellsInd,startCell);

            [sortedd, I] = sort(d);
            freeCellsInd = freeCellsInd(I);
  
            %mzaa verbundn sein 
            %stilll connected, also entweder ist mucilage mi weniger
            %konzentration oder direkter nachbar
            outerborder = sortedd < notConnectedEdgesValue * 1.1; +  %sortedd(1) + 
            %TODO hier vielleicht nochmal eukliduscher Abstand berechnne
            outerborderInd = freeCellsInd(outerborder);%heeeeeeree
            ind = outerborderInd;
            
            
end
function mucilageVector = updateMucilageVector()
indMucilagechanged = find( ( mucilageVector == 1 ) &  ( mucilageConcVector == 0 ) );
    if(~isempty(indMucilagechanged))
        mucilageVector(indMucilagechanged) = 0;
        bulkVector(indMucilagechanged) = 0;
        for i = 1:numel(indMucilagechanged)
            rootParticleList = rootParticleList(rootParticleList ~= indMucilagechanged(i));
        end
    end
end