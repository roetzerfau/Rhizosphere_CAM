function [porosity] = calculatePorosity(g,mantles,rootVector, soilMaterial)
%Reorganisation of rhizosphere soil pore structure by wild plant species in compacted soils 
%The emergent rhizosphere: imaging the development of the porous architecture at the root-soil interface
porosity = zeros(numel(mantles)-1,1);
N = g.NX;
geoSoil = reshape(soilMaterial, [N,N]);
geoRoot = reshape(rootVector, [N,N]);
D = bwdist(geoRoot);
for i = 1:numel(mantles)-1
    soilOI = geoSoil(D > mantles(i) & D <= mantles(i+1));
    porosity(i) = sum(1-soilOI,'all')/numel(soilOI);
end
end