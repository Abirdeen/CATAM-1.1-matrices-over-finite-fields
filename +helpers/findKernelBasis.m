function U = findKernelBasis(A, p, inverse)
    % FINDKERNELBASIS   Finds a basis for the kernel of a matrix.
    % 
    % parameters
    % ----------
    % A : Any matrix.
    % p : Prime modulus.
    % inverse : array of multiplicative inverses mod p, e.g. from FINDMODULARINVERSES.
    %
    % returns
    % -------
    % U : A basis of ker(A). The ith column is the ith basis element. If A has no kernel, U will be an n x 1 zero vector.
    %
    % examples
    % --------
    % U = FINDKERNELBASIS([1 2 0 3 0;0 0 1 4 0;0 0 0 0 1; 0 0 0 0 0], 5, [1,3,2,4]) returns [3 2; 1 0; 0 1; 0 1; 0 0]
    [~,n] = size(A);
    [B,~,pivots,rank] = helpers.gaussianElimination(A, p, inverse);
    if rank >= n
        U = zeros([n,1]);
        return
    end
    nonpivots = setdiff(1:n, pivots);
    U = zeros([n, n-rank]);
    for i = 1:numel(nonpivots)
        k = nonpivots(i);
        U(k,i) = 1;
        for s = 1:numel(pivots)
            j_s = pivots(s);
            if j_s > k
                break
            end
            U(j_s,i) = -B(s,k);
        end
    end
    U = mod(U,p);
end