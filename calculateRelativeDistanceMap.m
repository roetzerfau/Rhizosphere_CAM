function relDistanceMap = calculateRelativeDistanceMap(g,rootGraph,particleList, particleVector, rootInitialCellInd)
    
    %relDistanceMap = zeros(size(particleVector));
    borderpoints = getParticleSurface(g,particleList, particleVector);
    borderpoints = borderpoints{1};
    d_b = zeros(numel(particleVector),numel(borderpoints));
    for i = 1:numel(borderpoints)
        d = distances(rootGraph, borderpoints(i));
%         img = reshape(d, [g.NX g.NX]);
%         imshow(img, [])
        d_max = d(rootInitialCellInd);
        d = d./d_max;
        d(d>1) = 1;
        d_b(:,i) = d;
%         img = reshape(d, [g.NX g.NX]);
%         imshow(img, [])
    end
    if(numel(borderpoints) > 0)
        relDistanceMap = min(d_b,[],2);
    else
        relDistanceMap = zeros(size(particleVector));
    end
    relDistanceMap = 1-relDistanceMap;
%     img = reshape(relDistanceMap, [g.NX g.NX]);
%     imshow(img, [])
    
    
    %distance von jedem Borderpunkt
% dann geteilt durch distance zu mittelpunkt
%
%max wert 
end