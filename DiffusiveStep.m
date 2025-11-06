function [dfp,uDG,q1DG,q2DG] = DiffusiveStep(g,param,dfp,uDG,q1DG,q2DG,bulkVector)
zero        = @(x,y) x-x + 0;
    one         = @(x,y) x-x + 1;  
one_u         = @(x,y) x-x + 1 * param.DOC;
     if isempty(uDG)
            uDG         = projectAlg2DGsub(g, one_u, dfp.p, dfp.ord, dfp.hatMc, dfp.airVectorSub);
            q1DG        = projectAlg2DGsub(g, zero, dfp.p, dfp.ord, dfp.hatMc, dfp.airVectorSub);
            q2DG        = projectAlg2DGsub(g, zero, dfp.p, dfp.ord, dfp.hatMc, dfp.airVectorSub);
    else
        uDG         = uDG(dfp.airVectorSub,:);
        q1DG        = q1DG(dfp.airVectorSub,:);
        q2DG        = q2DG(dfp.airVectorSub,:);
    end

    sysU        = reshape(uDG', dfp.NT*dfp.N, 1);
    sysQ1       = reshape(q1DG', dfp.NT*dfp.N, 1);
    sysQ2       = reshape(q2DG', dfp.NT*dfp.N, 1);

    dfp.sysX = [    sysQ1                   ;   sysQ2                   ;   sysU                                ];
   
    
  
    uDG    = zeros(g.numTsub, dfp.N);
    uDG(dfp.airVectorSub,:) = reshape( dfp.sysX( 2 * dfp.NT * dfp.N + 1 : 3 * dfp.NT * dfp.N ), dfp.N, dfp.NT )';

    q1DG        = uDG;
    q2DG        = uDG;
  
   
    uDG(bulkVector == 1,1) = 0;

    uDG         = uDG;  



     if size(uDG,1) == g.numTsub
                uDG = uDG(dfp.airVectorSub,:);
     end
    sysU        = reshape(uDG', dfp.NT*dfp.N, 1);
    sysX = [    sysQ1                   ;   sysQ2                   ;   sysU                                ];
%         
    
    sysX = ( dfp.sysW + dfp.tau * dfp.sysA ) \ ( dfp.sysW * sysX + dfp.tau * dfp.sysV );

    
    uDG         = zeros(g.numTsub, dfp.N);
    uDG(dfp.airVectorSub,:) = reshape( sysX( 2 * dfp.NT * dfp.N + 1 : 3 * dfp.NT * dfp.N ), dfp.N, dfp.NT )';
% 
%         uDGbac      = zeros(g.numTsub, dfp.N);
%         uDGbac(dfp.airVectorSub,:) = reshape( sysXbac( 2 * dfp.NT * dfp.N + 1 : 3 * dfp.NT * dfp.N ), dfp.N, dfp.NT )';
% %        uDG(~dfp.airVectorSub,1) = -1;
%         
    %concAir = computeConcAir(g, uDG, dfp.ord);
end