function [ globG ] = assembleGlobGsub( g , hatGc , hatGx , hatGy , K11DG , ...
    K12DG , K21DG , K22DG , airVec , numAir )

[~, N] = size(K11DG);
K = numAir;
globG = cell(2,1);
globG{1} = sparse(K*N, K*N);
globG{2} = sparse(K*N, K*N);

for i = 1 : N
    globG{1} = globG{1} - kron( spdiags(K11DG(:,i) .* g.ACBDsub(airVec), 0, K, K), hatGx(:,:,i,1) ) ...
                        - kron( spdiags(K11DG(:,i) .* g.DAsub(airVec)  , 0, K, K), hatGc(:,:,i,1) ) ...
                        + kron( spdiags(K11DG(:,i) .* g.ACBDsub(airVec), 0, K, K), hatGy(:,:,i,2) ) ...
                        + kron( spdiags(K11DG(:,i) .* g.BAsub(airVec) - K21DG(:,i) .* g.deltaX, ...
                            0, K, K), hatGc(:,:,i,2) );
    globG{2} = globG{2} - kron( spdiags(K12DG(:,i) .* g.ACBDsub(airVec), 0, K, K), hatGx(:,:,i,1) ) ...
                        - kron( spdiags(K12DG(:,i) .* g.DAsub(airVec)  , 0, K, K), hatGc(:,:,i,1) ) ...
                        + kron( spdiags(K12DG(:,i) .* g.ACBDsub(airVec), 0, K, K), hatGy(:,:,i,2) ) ...
                        + kron( spdiags(K12DG(:,i) .* g.BAsub(airVec) - K22DG(:,i) .* g.deltaX, ...
                            0, K, K), hatGc(:,:,i,2) );
end  % for

end  % function