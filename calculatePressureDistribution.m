function pressureDistributionVector = calculatePressureDistribution(g, pressurePointsInd, bulkVector, noPressurePointsInd)

    pressureDistributionVector = zeros(g.numT, 1);
    pressurePointsConnectedBulk = zeros(g.numT, 1);
     pressurePointsConnectBulkInd = [];
    N = g.NX;
    bulkVector_bw = (reshape(bulkVector, [N N]));
    CC = bwconncomp(bulkVector_bw, 4);

   
    for p = 1:numel(pressurePointsInd)
		if(pressurePointsConnectedBulk(pressurePointsInd) == 1)
			pressurePointsConnectBulkInd = [pressurePointsConnectBulkInd pressurePointsInd];
			fprintf('isschon \n')
			continue;
		end
        oo = cellfun(@(m) m == pressurePointsInd(p),CC.PixelIdxList, 'UniformOutput', false);
        Match1 = cellfun(@sum, oo);
        Match = find(Match1 == 1);
        
        if(~isempty(Match))
            connectedBulkInd  =  CC.PixelIdxList{Match};    
            pressurePointsConnectedBulk(connectedBulkInd) = 1;
            pressurePointsConnectBulkInd = [pressurePointsConnectBulkInd pressurePointsInd];
        end
    end
     %-------------------------------------------------------
    % Pressure Distributoin - pressure Points ohne Verbindung keine
    % Auswirkung
    [X,Y] = meshgrid(1:N, 1:N);
    for i = 1:numel(pressurePointsConnectBulkInd)
        
        
        X_0 = X(pressurePointsConnectBulkInd(i));
        Y_0 = Y(pressurePointsConnectBulkInd(i));
        F = -((X-X_0).^2 + (Y-Y_0).^2) + (N/2)^2 ;
        F(F < 0) = 0;
        F = reshape(F, [N*N 1]);
        M = cat(2,pressureDistributionVector(:),F(:));
		pressureDistributionVector = reshape(max(M,[],2),size(pressureDistributionVector));
        %pressureDistributionVector = pressureDistributionVector + F;
    end
    pressureDistributionVector = normalize(pressureDistributionVector,'range',[0 1]);
    %---------------------------------------------------------------------------------
    % Combine Pressure mit Connected with Root
    pressureArea = pressurePointsConnectedBulk;
    pressureArea(pressurePointsInd) = 1;
    pressureArea(noPressurePointsInd) = 0;
    pressureDistributionVector = pressureDistributionVector.*pressureArea;
    %weiß nicht ob nötig:
    pressureDistributionVector(pressurePointsInd) = 1;
%     todo
%     pressureDistributionVector(pressureDistributionVector > 0) = normalize(pressureDistributionVector(pressureDistributionVector > 0),'range',[0 1]);
    
    
    
    
end
