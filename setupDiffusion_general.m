function dfp = setupDiffusion_general(g)

dfp.tau         = 1;                        
dfp.epsilon     = 1e-6;                     % Everything < epsilon is defined as zero

dfp.p           = 1;                        % Order of Polynomials for LDG Disc.
dfp.ord         = 4;                        % Order of Gaussian Integration Rule
dfp.eta         = 1;                        % Penalty Parameter for LDG

computeBasesOnQuad(dfp.p, dfp.ord);

[dfp.hatMc, dfp.hatMx]          = computeHatM(dfp.p, dfp.ord);
[dfp.hatGc, dfp.hatGx, dfp.hatGy]   = computeHatG(dfp.p, dfp.ord);
[dfp.hatHc, dfp.hatHx, dfp.hatHy]   = computeHatH(dfp.p, dfp.ord);
dfp.hatRdiag                = computeHatRdiag(dfp.p, dfp.ord);
dfp.hatRoffdiag             = computeHatRoffdiag(dfp.p, dfp.ord);
dfp.hatSdiag                = computeHatSdiag(dfp.p, dfp.ord);
dfp.hatSoffdiag             = computeHatSoffdiag(dfp.p, dfp.ord);

end