function concentrationVector = spreadConcentration(g, concentrationVector, restrictVector, maxValueVector, minValueVector)
previousConcentration = sum(concentrationVector);  
exceedMaximium_ind = find(concentrationVector > (maxValueVector + minValueVector));
    
    for i = 1:numel(exceedMaximium_ind)
        overshoot_candidate = exceedMaximium_ind(i);
        
        isOverlapped = restrictVector(overshoot_candidate) > 0;
        if(isOverlapped)
            overshoot_value = concentrationVector(overshoot_candidate);
        else
            overshoot_value =  concentrationVector(overshoot_candidate)/2;% - maxValueVector(overshoot_candidate);
        end
        
        %candidates = candidates(randperm(length(candidates)));
        candidates = overshoot_candidate;
        layer = 1;
        previous_candidates = overshoot_candidate;
        while(overshoot_value > 0)    
            
             if(isOverlapped)
                candidates = stencil( g.NX , g.NX , candidates , layer );
             else
                candidates = unique(stencil( g.NX , g.NX , candidates , 1 ));
                 
                isOccupied = restrictVector(candidates) > 0;
                candidates(isOccupied) = [];
            end
            
            
            
            
            undershoot_candidates = candidates(concentrationVector(candidates) < maxValueVector(candidates)); 
            undershoot_candidates = undershoot_candidates(randperm(length(undershoot_candidates)));
            
            for c = reshape(undershoot_candidates,1,[])
                capacity = maxValueVector(c) - concentrationVector(c);
                fraction_2_spread = min([capacity, overshoot_value]);
                concentrationVector(c)= concentrationVector(c) + fraction_2_spread;
                concentrationVector(overshoot_candidate) = concentrationVector(overshoot_candidate) - fraction_2_spread;
                overshoot_value = overshoot_value -  fraction_2_spread;
            end
            
           if(isequal(previous_candidates,candidates))
               equalConcentration = sum(concentrationVector(candidates))/numel(candidates);
               concentrationVector(candidates) = equalConcentration;
               break;
           end
            layer = layer +1;
            previous_candidates = candidates;
        end
        
        
   
    end
    currentConcentration = sum(concentrationVector);  
    if(abs(previousConcentration - currentConcentration) > 0.000001)
              fprintf("falsch spread")
              %error('Falsch spread', abs(previousConcentration - currentConcentration))
    end
end
%Wenn mehr und es gibt keinen freien platz mehr -> gleich verteilt
%biomasse darf nicht wachsen wenn es maximale capacity überschreitet