function [particleSurfaceList] = getParticleSurface(g, particleList, particleTypeVector)
    particleSurfaceList = cell(length(particleList),1);
    for i = 1 : length(particleList)
           for j = 1 : length(particleList{i})
               ind = particleList{i}(j);
               possibleNeighbors = stencil( g.NX, g.NX, ind, 1); 
               possibleNeighbors = possibleNeighbors(2:end);
               for neigh = 1 : 4
                   if (particleTypeVector(possibleNeighbors(neigh)) ~= particleTypeVector(ind) )
                        particleSurfaceList{i} = [particleSurfaceList{i} ind];
                        break;
                   end
               end
           end
    end
end