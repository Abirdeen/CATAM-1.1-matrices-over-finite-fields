function [B, R, rank] = gaussianElimination(A, p, inverse)
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
    %
    % examples
    % --------
    % [B, R, rank] = GAUSSIANELIMINATION([1,2;1,1;1,0], 3, []) returns [1 4 5 2 3 6]
    % If a number has no inverse, returns 0
    % INVERSES = FINDMODULARINVERSES(4) returns [1 0 3]
    [n,m] = size(A);
    B = mod(A,p);
    R = eye(n);
    row_transps = 0;
    for i = 1:m
        nonzero_index = 0;
        % Search for a non-zero index in column i
        for j = row_transps+1:n
            if B(j,i) == 0
                continue
            end
            nonzero_index = j;
            break
        end
        % If no non-zero index is found, move on to next column
        if nonzero_index == 0
            continue
        end
        % Move the non-zero row up
        % Normalise the row (so the first non-zero index is 1)
        row_transps = row_transps + 1;
        Q = T(n,row_transps,nonzero_index);
        Q = D(n,row_transps,B(nonzero_index,i), inverse)*Q;
        B = mod(Q*B,p);
        R = mod(Q*R,p);
        % Eliminate the column from other rows
        for j = 1:n 
            if j == row_transps || B(j,i) == 0
                continue
            end
            Q = S(n, row_transps, B(j,i), j);
            B = mod(Q*B, p);
            R = mod(Q*R, p);
        end
    end
    rank = row_transps;
end

function R = T(n,i,j)
    R = eye(n);
    if i>n || j>n
        return 
    end
    [R(i,i),R(j,j), R(i,j),R(j,i)] = deal(0,0,1,1);
end

function R = D(n, i, a, inverse)
    R = eye(n);
    if i>n
        return
    end
    R(i,i) = inverse(a);
end

function R = S(n, i, a, j)
    R = eye(n);
    if i>n || j>n
        return
    end
    R(j,i) = -a;
end