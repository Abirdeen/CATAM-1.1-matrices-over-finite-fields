function A_inverse = MoorePenroseInverse(A, p, inverse)
    % INVERTMATRIX   Find A^-1 mod p.
    % 
    % parameters
    % ----------
    % A : A square matrix to invert.
    % p : A prime number.
    % inverse : Array of multiplicative inverses mod p, e.g. from FINDMODULARINVERSES.
    %
    % returns
    % -------
    % A_inverse : A^+, the Moore-Penrose pseudo-inverse of A. If A is square, this is the usual inverse. If 
    % A does not have full rank, returns a matrix of zeros with the same dimensions as (the transpose of) A.
    %
    % examples
    % --------
    % A_inverse = INVERTMATRIX([0,1;2,0],3,[1,2]) returns [0,2;1,0]
    % A_inverse = INVERTMATRIX([0,1,0;4,1,1],5,[1,3,2,4]) returns [3,2;1,0;2,3]
    arguments
          A {mustBeMatrix}
          p {mustBePositive}
          inverse {mustBeRow}
    end
    [m,n] = size(A);
    if ~(m==n)
          B = mod(A*transpose(A),p);
    else
          B = mod(A,p);
    end
    det_B = mod(det(B),p);
    if det_B==0
          A_inverse = zeros([n,m]);
          return
    end
    B_inverse = mod(inverse(det_B).*adjoint(B),p);
    if ~(m==n)
          A_inverse = mod(transpose(A)*B_inverse,p);
    else
          A_inverse = B_inverse;
    end
end