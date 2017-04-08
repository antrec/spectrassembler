# Julia script to compute eigenvector of Laplacian

# Load indexes and construct sparse matrix
iis = readcsv(open(ARGS[1]), Int32)
jjs = readcsv(open(ARGS[2]), Int32)
vvs = readcsv(open(ARGS[3]), Int32)
A = sparse(iis[:]+1, jjs[:]+1, vvs[:])

# Compute Laplacian
Lap = spdiagm(sum(A,2)[:,1], 0, size(A)[1], size(A)[2]) - A + 1e-9*speye(Float64, size(A)[1]);

# Get largest eigenvalue
(lbdamax, ~) = eigs(Lap, nev=1);

# Construct lambda_max*I - Laplacian, its second largest eigenvalue is the fiedler value we are looking for
negLap = lbdamax[1]*speye(Float64, size(Lap)[1]);
negLap = negLap - Lap;

# Get eigenvalue and eigenvector
(lbdadiff, V2) = eigs(negLap, nev=6, ncv=20, tol=1.e-9, v0=ones((size(Lap)[1],)), maxiter=10000000);
fidval = lbdamax[1] - lbdadiff[2]
fidvec = V2[:,2];

# Sort to recover permutation (spectral ordering)
permu = sortperm(fidvec)

writecsv(ARGS[4], permu')
