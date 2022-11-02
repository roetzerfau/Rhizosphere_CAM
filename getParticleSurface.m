function [particleSurfaceList] = getParticleSurface(g, particleList, particleTypeVector)
    particleSurfaceList = [];
   for j = 1 : length(particleList)
       ind = particleList(j);
       possibleNeighbors = stencil( g.NX, g.NX, ind, 1); 
       possibleNeighbors = possibleNeighbors(2:end);
       for neigh = 1 : 4
           if (particleTypeVector(possibleNeighbors(neigh)) ~= particleTypeVector(ind) )
                particleSurfaceList = [particleSurfaceList ind];
                break;
           end
       end
   end
end
