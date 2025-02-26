function [B, R, pivots, rank] = gaussianElimination(A, p, inverse)
    % GAUSSIANELIMINATION   Finds the echelon-reduced form of A.
    % 
    % parameters
    % ----------
    % A : Matrix to reduce.
    % p : Prime modulus.
    % inverse : array of multiplicative inverses mod p, e.g. from FINDMODULARINVERSES.
    %
    % returns
    % -------
    % B : Echelon-reduced form of A.
    % R : Matrix such that B = R*A.
    % rank : The rank of A.
    % pivots : Array of column indexes of pivots of B.
    %
    % examples
    % --------
    % [B, R, pivots, rank] = HELPERS.GAUSSIANELIMINATION([1,2;1,1;1,0], 3, [1,2]) returns
    % B = 1 0       R = 2 2 0       pivots = 1      rank = 2
    %     0 1           1 2 0                2
    %     0 0           1 1 1
    arguments
        A {mustBeMatrix}
        p {mustBePositive}
        inverse {mustBeRow}
    end
    [m,n] = size(A);
    B = mod(A,p);
    R = eye(m);
    row_transps = 0;
    pivots = zeros(m);
    for j = 1:n
        nonzero_index = 0;
        % Search for a non-zero index in column i
        for i = row_transps+1:m
            if B(i,j) == 0
                continue
            end
            nonzero_index = i;
            break
        end
        % If no non-zero index is found, move on to next column
        if nonzero_index == 0
            continue
        end
        % Move the non-zero row up
        % Normalise the row (so the first non-zero index is 1)
        row_transps = row_transps + 1;
        Q = T(m,row_transps,nonzero_index);
        Q = D(m,row_transps,B(nonzero_index,j), inverse)*Q;
        B = mod(Q*B,p);
        R = mod(Q*R,p);
        % Add column to pivots
        pivots(row_transps) = j;
        % Eliminate the column from other rows
        for i = 1:m 
            if i == row_transps || B(i,j) == 0
                continue
            end
            Q = S(m, row_transps, B(i,j), i);
            B = mod(Q*B, p);
            R = mod(Q*R, p);
        end
    end
    rank = row_transps;
    pivots = nonzeros(pivots);
end

function R = T(m,j,k)
    R = eye(m);
    if j>m || k>m
        return 
    end
    [R(j,j),R(k,k), R(j,k),R(k,j)] = deal(0,0,1,1);
end

function R = D(m, j, a, inverse)
    R = eye(m);
    if j>m
        return
    end
    R(j,j) = inverse(a);
end

function R = S(m, j, a, k)
    R = eye(m);
    if j>m || k>m
        return
    end
    R(k,j) = -a;
end