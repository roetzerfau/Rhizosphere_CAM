function distance_map = distanceMap(bulkVector, referenceCell)
    distance_map = bulkVector;
    distance_map(:) = -1;
    NX = sqrt(length(bulkVector));
    waitList = [];
    waitList = [waitList, referenceCell];
    
    distance = 0;
    distance_map(referenceCell)= distance;
    while isempty(waitList) == 0
        distance = distance +1;
        currentBulk = waitList(1);
        waitList = waitList(2:end);
        %if bulkVectorFlag(currentBulk) < 0
        %    bulkVectorFlag(currentBulk) = particleCounter;
        %end            
        possibleNeighbors = stencil( NX, NX, currentBulk, 1); 
        possibleNeighbors = possibleNeighbors(2:end);
        neighbors = possibleNeighbors(bulkVector(possibleNeighbors) == 1);
        for neighborBulk = neighbors

            if distance_map(neighborBulk) <= -1
            distance_map(neighborBulk) = distance_map(currentBulk) + 1;
            waitList = [waitList, neighborBulk];
            end
        end
    end



end