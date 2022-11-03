function [relDistanceMap, borderpointsVector]  = calculateRelativeDistanceMap(g,graph,particleList, particleVector, referenceCellInd)
    
    %relDistanceMap = zeros(size(particleVector));
    borderpoints = getParticleSurface(g,particleList, particleVector);
    borderpointsVector = zeros(g.numT, 1);
    borderpointsVector(borderpoints) = 1;
    d_b = zeros(numel(particleVector),numel(borderpoints));
    for i = 1:numel(borderpoints)
        d = distances(graph, borderpoints(i));
%         img = reshape(d, [g.NX g.NX]);
%         imshow(img, [])
        d_max = d(referenceCellInd);
        d = d./d_max;
        d(d>1) = Inf;
        d_b(:,i) = d;
%         img = reshape(d, [g.NX g.NX]);
%         imshow(img, [])
    end
    if(numel(borderpoints) > 0)
        relDistanceMap = min(d_b,[],2);
    else
        relDistanceMap = zeros(size(particleVector));
    end
    relDistanceMap(~particleVector) = Inf;
    relDistanceMap = 1-relDistanceMap;
    %Problem Periodizität Lösung:
    
    
    
    %relDistanceMap(borderpoints) = 1;
%     img = reshape(relDistanceMap, [g.NX g.NX]);
%     imshow(img, [])
    
    
    %distance von jedem Borderpunkt
% dann geteilt durch distance zu mittelpunkt
%
%max wert 
end