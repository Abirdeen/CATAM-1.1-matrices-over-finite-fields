clc 
clear 
close all

primes = [2 3 5 7 11 13 17 19 23 29 31 37 41 43 47 48 53 59 61 67 71 73 79 83 89 97];

A = [11 1 7 2 0
     8  0 2 5 11
     2  1 2 6 5
     7  4 5 3 1 ];

B = [0 1 1 3 5  2
     1 2 3 8 9  0
     0 1 1 2 3  2
     2 1 3 7 9  1
     2 1 3 8 10 0];

%C = [1,1;1,2];

%D = [1 2 0 3 0; 0 0 1 4 0; 0 0 0 0 1; 0 0 0 0 0];

%p = 7;
%inverse = helpers.findModularInverses(p);

for i = 1:10
    p = primes(i);
    disp(p)
    inverse = helpers.findModularInverses(p);

    U = helpers.findKernelBasis(B, p, inverse);
    disp(U)
end

