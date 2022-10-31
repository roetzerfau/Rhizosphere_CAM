function [mucilageConcVector, mucilageVector, mucilageParticleList] = ...
    calculateMucilageDecay(g, parameters, bulkVector, mucilageVector, mucilageParticleList, mucilageConcVector)
    
    tau = 1;
    for i = 1 : length(mucilageParticleList)
        for j = 1 : length(mucilageParticleList{i})
            ind = mucilageParticleList{i}(j);

            st = stencil( g.NX, g.NX, ind(i), 1);

            concOld = mucilageConcVector(ind);
            l = sum(bulkVector(st(2:end)));
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
   end
    
%     indMucilagechanged = find( ( mucilageVector == 1 ) &  ( mucilageConcVector == 0 ) );
%     
%     mucilageVector(indMucilagechanged) = 0;
%     mucilageParticleList = mucilageParticleList(mucilageParticleList ~= indMucilagechanged);
%growing rate verkleiner
%homgenisieren durch mittelwert aus nachbarn
end