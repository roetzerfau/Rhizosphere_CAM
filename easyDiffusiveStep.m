function concentrationVector = easyDiffusiveStep(g, concentrationVector, restrictVector, range)

    previousConcentration = sum(concentrationVector);
    Concentration_Phase = find(concentrationVector > 0);
    Concentration_Phase = Concentration_Phase(randperm(length(Concentration_Phase)));
    
    Concentration_Phase_Occupied = Concentration_Phase(restrictVector(Concentration_Phase) > 0 );
    Concentration_Phase_NotOccupied = Concentration_Phase(restrictVector(Concentration_Phase) == 0);
    
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
% concentrationVector(Concentration_Phase_Occupied(i))
% numel(diffusive_area)
        equalConcentration = concentrationVector(Concentration_Phase_Occupied(i))/numel(diffusive_area);
        concentrationVector(Concentration_Phase_Occupied(i)) = 0;
        concentrationVector(diffusive_area) =concentrationVector(diffusive_area)+ equalConcentration;
    end
    
    for i = 1:numel(Concentration_Phase_NotOccupied)
        diffusive_area = Concentration_Phase_NotOccupied(i);
        
        for r = 1:range
            diffusive_area = unique(stencil( g.NX , g.NX , diffusive_area , 1 ));
            isOccupied = restrictVector(diffusive_area) > 0;         
            diffusive_area(isOccupied) = [];
        end
        if(numel(diffusive_area) == 0)
            printf('check')
  
        end
        equalConcentration = sum(concentrationVector(diffusive_area))/numel(diffusive_area);
        concentrationVector(diffusive_area) = equalConcentration;
    end
    
    
    

%     for i = 1:numel(Concentration_Phase_NotOccupied)
%         %concentration_area = Concentration_Phase(i);  
%         concentration_area = unique(stencil( g.NX , g.NX , Concentration_Phase_NotOccupied(i) , 5 ));
%         diffusive_area = Concentration_Phase_NotOccupied(i);
%         
%         for r = 1:range
%             diffusive_area = unique(stencil( g.NX , g.NX , diffusive_area , 1 ));
%             isOccupied = restrictVector(diffusive_area) > 0;         
%             diffusive_area(isOccupied) = [];
%         end
%         if(numel(diffusive_area) == 0)
%             printf('check')
%             
%         end
%         equalConcentration = sum(concentrationVector(concentration_area))/numel(diffusive_area);
%         concentrationVector(diffusive_area) = equalConcentration;
%     end
    currentConcentration = sum(concentrationVector);  
    if(abs(previousConcentration - currentConcentration) > 0.0001)
              abs(previousConcentration - currentConcentration)
              error('Falsch diffusive', abs(previousConcentration - currentConcentration))
    end
    if(sum(restrictVector(concentrationVector > 0)) ~= 0)
         sum(restrictVector(concentrationVector > 0))
         error('Falsch diffusive occupied')
    end
end