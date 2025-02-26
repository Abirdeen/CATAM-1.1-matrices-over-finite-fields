function [U,W,U_plus_W,U_cap_W] = secondIsomorphismTheoremBases(A, B, p, inverse)
    % SECONDISOMORPHISMTHEOREMBASES   Finds bases for U, W, U+W and UnW
    % 
    % parameters
    % ----------
    % A, B : Matrices with row spaces U and W respectively.
    % p : Prime modulus.
    % inverse : array of multiplicative inverses mod p, e.g. from FINDMODULARINVERSES.
    %
    % returns
    % -------
    % U, W, U_plus_W, U_cap_W : Bases of U, W, U+W and UnW as columns of matrices.
    %
    % examples
    % --------
    % [U, W, U_plus_W, U_cap_W] = HELPERS.SECONDISOMORPHISMTHEOREMBASES([1,2,1;1,1,1], [0,0,1], 3, [1,2]) returns
    % U = 1 0       W = 0       U_plus_W = 1 0 0      U_cap_W = 0
    %     0 1           0                  0 1 0                0
    %     1 0           1                  0 0 1                0
    arguments
        A {mustBeMatrix}
        B {mustBeMatrix}
        p {mustBePositive}
        inverse {mustBeRow}
    end
    [AA,~,~, rank_A] = helpers.gaussianElimination(A,p,inverse);
    [BB,~,~, rank_B] = helpers.gaussianElimination(B,p,inverse);
    U = transpose(AA(1:rank_A, :));
    W = transpose(BB(1:rank_B, :));

    [AB,~,~,rank_AB] = helpers.gaussianElimination([A;B],p,inverse);
    U_plus_W = transpose(AB(1:rank_AB, :));

    ann_U = helpers.findKernelBasis(A,p,inverse);
    ann_W = helpers.findKernelBasis(B,p,inverse);
    U_cap_W = helpers.findKernelBasis([transpose(ann_U); transpose(ann_W)],p,inverse);
end