function [dfp] = setupDiffusion_specific_domain(g,dfp, bulkVector)

    zero        = @(x,y) x-x + 0;
    one         = @(x,y) x-x + 1;
    one_diff    = @(x,y) x-x + 1 * 1.94 * 10^-10 * 10^8;%6.73 * 10^-6 * 10^8; 1.94 * 10^-10 * 10^8;
  
    markEbdr        = g.idE == 8 | g.idE == 8 | g.idE == 8 | g.idE == 8 | g.idE == 8;

    idE                  = zeros(g.numE, 1);
    idE(g.E0T(bulkVector == 1,1)) = idE(g.E0T(bulkVector == 1,1)) + 1;
    idE(g.E0T(bulkVector == 1,2)) = idE(g.E0T(bulkVector == 1,2)) + 1;
    idE(g.E0T(bulkVector == 1,3)) = idE(g.E0T(bulkVector == 1,3)) + 1;
    idE(g.E0T(bulkVector == 1,4)) = idE(g.E0T(bulkVector == 1,4)) + 1;
    idE = idE + markEbdr;
    
    idEneum         = idE == 1;
    idEairy         = idE == 0;
   
    markE0Tneum     = idEneum(g.E0T);
    markE0Tairy     = idEairy(g.E0T);
    airVector   = ~bulkVector;
    
    markE0TneumSub = markE0Tneum(1:g.numTsub,:);
    markE0TairySub = markE0Tairy(1:g.numTsub,:);
    dfp.airVectorSub = airVector(1:g.numTsub);
    markE0TneumSub = markE0TneumSub(dfp.airVectorSub,:);
    markE0TairySub = markE0TairySub(dfp.airVectorSub,:);
    numAir = norm(dfp.airVectorSub.*ones(size(dfp.airVectorSub)),1);

    %% Assembling global Matrices (independedfp.NT of u, q, concedfp.NTrations)
    globM       = assembleGlobMsub(g, dfp.hatMc, dfp.hatMx, dfp.airVectorSub, numAir);
    globH       = assembleGlobHsub(g, dfp.hatHc, dfp.hatHx, dfp.hatHy, dfp.airVectorSub, numAir);
    globQ       = assembleGlobQsub(g, markE0TairySub, dfp.hatSdiag, dfp.hatSoffdiag, dfp.airVectorSub, numAir);
    globQN      = assembleGlobQNsub(g, markE0TneumSub, dfp.hatSdiag, dfp.airVectorSub, numAir);
    globS       = assembleGlobSsub(g, markE0TairySub, dfp.hatSdiag, dfp.hatSoffdiag, dfp.eta, dfp.airVectorSub, numAir);


    fDG         = projectAlg2DGsub(g, zero, dfp.p, dfp.ord, dfp.hatMc, dfp.airVectorSub);
    K11DG       = projectAlg2DGsub(g, one_diff, dfp.p, dfp.ord, dfp.hatMc, dfp.airVectorSub);
    K12DG       = projectAlg2DGsub(g, zero, dfp.p, dfp.ord, dfp.hatMc, dfp.airVectorSub);
    K21DG       = projectAlg2DGsub(g, zero, dfp.p, dfp.ord, dfp.hatMc, dfp.airVectorSub);
    K22DG       = projectAlg2DGsub(g, one_diff, dfp.p, dfp.ord, dfp.hatMc, dfp.airVectorSub);

    dfp.NT = numAir;
    dfp.N = (dfp.p+1)^2;

   
   
    
    globL       = globM * reshape(fDG', dfp.NT*dfp.N, 1);
    globG       = assembleGlobGsub(g, dfp.hatGc, dfp.hatGx, dfp.hatGy, K11DG, K12DG, K21DG, K22DG, dfp.airVectorSub, numAir);
    globR       = assembleGlobRsub(g, markE0TairySub, dfp.hatRdiag, dfp.hatRoffdiag, K11DG, K12DG, K21DG, K22DG, dfp.airVectorSub, numAir);
    globKN      = assembleGlobKNsub(g, dfp.p, dfp.ord, markE0TneumSub, zero, dfp.airVectorSub, numAir);
  

   

    dfp.sysW = [    sparse(dfp.NT*dfp.N,dfp.NT*dfp.N)       ,   sparse(dfp.NT*dfp.N,dfp.NT*dfp.N)       ,   sparse(dfp.NT*dfp.N,dfp.NT*dfp.N)                   ;
                sparse(dfp.NT*dfp.N,dfp.NT*dfp.N)       ,   sparse(dfp.NT*dfp.N,dfp.NT*dfp.N)       ,   sparse(dfp.NT*dfp.N,dfp.NT*dfp.N)                   ;
                sparse(dfp.NT*dfp.N,dfp.NT*dfp.N)       ,   sparse(dfp.NT*dfp.N,dfp.NT*dfp.N)       ,   globM                               ];
    dfp.sysA = [    globM                   ,   sparse(dfp.NT*dfp.N,dfp.NT*dfp.N)       ,   globH{1} + globQ{1} + globQN{1}     ;
                sparse(dfp.NT*dfp.N,dfp.NT*dfp.N)       ,   globM                   ,   globH{2} + globQ{2} + globQN{2}     ;
                globG{1} + globR{1}     ,   globG{2} + globR{2}     ,   globS                               ];

    dfp.sysV = [    zeros(size(globM,1),1)  ;   zeros(size(globM,1),1)  ;   globKN + globL                      ];




   
   


end