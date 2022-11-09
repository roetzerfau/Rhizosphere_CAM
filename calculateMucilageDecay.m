function [mucilageConcVector, mucilageVector, rootParticleList, bulkVector] = ...
    calculateMucilageDecay(g, parameters, bulkVector, mucilageVector, rootParticleList, mucilageConcVector)
    
    tau = 1;
    mucilageParticleList = rootParticleList(mucilageVector(rootParticleList) == 1);
    for j = 1 : length(mucilageParticleList)
        ind = mucilageParticleList(j);

        st = stencil( g.NX, g.NX, ind, 1);

        concOld = mucilageConcVector(ind);
        decayRate = parameters.mucilageDecayRate * (4-sum(bulkVector(st(2:end))))/4;
        particleVector = bulkVector-mucilageVector;
        if( sum(particleVector(st(2:end))) >= 1)
            decayRate = decayRate * 0.3;
        end

        concNew = concOld *  exp(- decayRate * tau);

        mucilageConcVector(ind) = concNew;

        if(mucilageConcVector(ind) < parameters.minConcMucilage)
            mucilageConcVector(ind) = 0;
        end
    end
  
    
    indMucilagechanged = find( ( mucilageVector == 1 ) &  ( mucilageConcVector == 0 ) );
    if(~isempty(indMucilagechanged))
        mucilageVector(indMucilagechanged) = 0;
        bulkVector(indMucilagechanged) = 0;
        for i = 1:numel(indMucilagechanged)
            rootParticleList = rootParticleList(rootParticleList ~= indMucilagechanged(i));
        end
    end
%growing rate verkleiner
%homgenisieren durch mittelwert aus nachbarn
end