function [mucilageConcVector] = ...
    calculateMucilageDecay(g, mucilageDecayRate, bulkVector, mucilageVector, mucilageConcVector)
    
    tau = 1;
    mucilageParticleList = find(mucilageVector == 1);
    for j = 1 : length(mucilageParticleList)
        ind = mucilageParticleList(j);

        st = stencil( g.NX, g.NX, ind, 1);

        concOld = mucilageConcVector(ind);
        neigh = (bulkVector + mucilageVector) >= 1;
        decayRate = mucilageDecayRate * (4-sum( neigh(st(2:end))))/4;
        particleVector = bulkVector;
        if( sum(particleVector(st(2:end))) >= 1)
            decayRate = decayRate * 0.3;
        end

        concNew = concOld *  exp(- decayRate * tau);

        mucilageConcVector(ind) = concNew;
    end
  
    
%     indMucilagechanged = find( ( mucilageVector == 1 ) &  ( mucilageConcVector == 0 ) );
%     if(~isempty(indMucilagechanged))
%         mucilageVector(indMucilagechanged) = 0;
%         bulkVector(indMucilagechanged) = 0;
%         for i = 1:numel(indMucilagechanged)
%             rootParticleList = rootParticleList(rootParticleList ~= indMucilagechanged(i));
%         end
%     end
%growing rate verkleiner
%homgenisieren durch mittelwert aus nachbarn
end