function concentrationVector = easyDiffusiveStep_occupied(g, concentrationVector, restrictVector)

    previousConcentration = sum(concentrationVector);
    Concentration_Phase = find(concentrationVector > 0);
    Concentration_Phase = Concentration_Phase(randperm(length(Concentration_Phase)));
    
    Concentration_Phase_Occupied = Concentration_Phase(restrictVector(Concentration_Phase) > 0);
    Concentration_Phase_Occupied = [Concentration_Phase_Occupied; find(concentrationVector > 0.01)];
   
  

    for i = 1:numel(Concentration_Phase_Occupied)
        diffusive_area = [];
        layer = 1;
        distribute = false;
        while(~distribute)
            diffusive_area = unique(stencil( g.NX , g.NX , Concentration_Phase_Occupied(i) , layer ));
            isOccupied = restrictVector(diffusive_area) > 0;         
            diffusive_area(isOccupied) = [];
            layer = layer +1;
            if(numel(diffusive_area) ~= 0)
                distribute = true;
            end
        end

        equalConcentration = concentrationVector(Concentration_Phase_Occupied(i))/numel(diffusive_area);
        concentrationVector(Concentration_Phase_Occupied(i)) = 0;
        concentrationVector(diffusive_area) =concentrationVector(diffusive_area)+ equalConcentration;
    end
    
 
currentConcentration = sum(concentrationVector);  
    if(abs(previousConcentration - currentConcentration) > 0.000001)
              abs(previousConcentration - currentConcentration)
              error('Falsch diffusive %f \n', abs(previousConcentration - currentConcentration))
    end
    if(sum(restrictVector(concentrationVector > 0)) ~= 0)
         sum(restrictVector(concentrationVector > 0))
         error('Falsch diffusive occupied')
    end



end