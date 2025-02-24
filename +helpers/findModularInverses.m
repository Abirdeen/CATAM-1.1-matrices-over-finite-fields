function inverse = findModularInverses(p)
    % FINDMODULARINVERSES   Find inverses to integers mod p.
    % 
    % parameters
    % ----------
    % p : A prime number.
    %
    % returns
    % -------
    % inverse : an array of integers such that mod(inverse(a)*a,p) = 1.
    %
    % examples
    % --------
    % INVERSES = FINDMODULARINVERSES(7) returns [1 4 5 2 3 6]
    % If a number has no inverse, returns 0
    % INVERSES = FINDMODULARINVERSES(4) returns [1 0 3]
    arguments
        p {mustBePositive}
    end
    inverse = zeros(1,p-1);
    inverse(1) = 1;
    inverse(p-1) = p-1;
    for i = 2:p-2
        if not (inverse(i) == 0)
            continue
        end
        for j = i:p-1
            if mod(i*j,p) == 1
                inverse(i) = j;
                inverse(j) = i;
            end
        end
    end
end