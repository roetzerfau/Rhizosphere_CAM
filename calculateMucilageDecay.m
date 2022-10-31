function [mucilageConcVector, mucilageVector, mucilageParticleList] = ...
    calculateMucilageDecay(g, paramters, bulkVector, mucilageVector, mucilageParticleList, mucilageConcVector)
    
    tau = 1;
    for i = 1:numel(mucilageParticleList)

        st = stencil( g.NX, g.NX, mucilageParticleList(i), 1);
        
        concOld = mucilageConcVector(mucilageParticleList(i));
        decayRate = parameters.mucilageDecayRate * (4-sum(bulkVector(st(2:end))))/4;
        concNew = concOld *  exp(- decayRate * tau);
        
        mucilageConcVector(mucilageParticleList(i)) = concNew;
        
        if(mucilageConcVector(mucilageParticleList(i)) < paramters.minConcMucilage)
            mucilageConcVector(mucilageParticleList(i)) = 0;
        end
    end
    
    indMucilagechanged = find( ( mucilageVector == 1 ) &  ( mucilageConcVector == 0 ) );
    
    mucilageVector(indMucilagechanged) = 0;
    mucilageParticleList = mucilageParticleList(mucilageParticleList ~= indMucilagechanged);
%growing rate verkleiner
%homgenisieren durch mittelwert aus nachbarn
end