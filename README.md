# CATAM-1.1-matrices-over-finite-fields

A rewrite of my Cambridge CATAM undergraduate project, "Matrices over finite fields".

## MATLAB

This project is written in [MATLAB](https://www.mathworks.com/products/matlab.html), a proprietary numeric computing language, originally designed for working with matrices. The remarkably light syntax and fast matrix manipulation make it well-suited for this matrix-oriented project.

## The project

The aim of this project is to study algebraic invariants of vector spaces over the [finite fields](https://en.wikipedia.org/wiki/Finite_field) $GF(p)$, where $p$ is some prime number. It is well-known that $GF(p) \cong \mathbb{Z}_p$, the integers modulo $p$.

Our examples use small values of $p$, i.e., $p < 30$, and small matrixes, i.e. at most $10 \times 10$, but the methods are theoretically applicable in full generality.

### Modular inverses

It follows from [Bézout's lemma](https://en.wikipedia.org/wiki/B%C3%A9zout%27s_identity) that, for any $a$ and $n$ coprime, there exists a $b$ such that $ab \equiv 1 \mod n$. In particular, when working modulo a prime $p$, every non-zero number has a multiplicative inverse.

There are multiple algorithms to find the multiplicative inverse of $a$ modulo $p$. For computing a single inverse, the [Extended Euclid's algorithm](https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm) is an efficient option, with runtime $O(\ln a )$. One might also use the [Fermat-Euler Theorem](https://en.wikipedia.org/wiki/Euler%27s_theorem), and compute $a^{-1}=a^{p-2}$.

We find it useful to compute all inverses simultaneously. Doing this via the Extended Euclid's algorithm would have a runtime of $O(\ln 1 + \ln 2 + ... + \ln (p-1)) = O(\ln(p!))=O(p \ln p)$. A more naive algorithm - checking if $ab \equiv 1 \mod p$ for every $1 < a < b < p-1$ - has a slower runtime of $O(p^2)$, but is easier to implement. Since we work with only small values of $p$, we use the latter method, with some minor optimisations.

### Reduced row echelon form

Every matrix can be put in [reduced row echelon form](https://en.wikipedia.org/wiki/Row_echelon_form) by a sequence of elementary operations - swapping rows, multiplying a row by a constant, or subtracting a multiple of one row from another. The process of putting a matrix in this form is called [Gaussian elimination](https://en.wikipedia.org/wiki/Gaussian_elimination).

From the reduced row echelon form, one can read off the [rank](https://en.wikipedia.org/wiki/Rank_(linear_algebra)), as well as a basis of the row space - since the elementary operations are automorpisms of the row space, a basis is just the non-zero rows, and the rank is the number of such rows.

## Problems

The original CATAM project involved certain explicit questions and problems, which are reproduced (and solved) here.

### Problem 1:

Write a program to store the inverses of the non-zero elements of $\mathbb{Z}_p$ in an array of length $p − 1$. Find the inverses by testing, for each $a$ in the range $1 < a < p − 1$, all values of $b$ in the range $1 < b < p − 1$ until you find one which works and then store it. Describe any very simple modification to speed up this procedure (say by a factor of 2).

#### Solution:

This program is implemented as `helpers.findModularInverses`. A number of small efficiencies are implemented:

- It is always true that $1^2 \equiv (p-1)^2 \equiv 1 \mod p$, so these $1^{-1} = 1$ and $(p-1)^{-1} = p-1$ are stored at the start. These are the only values for which $a^2 \equiv 1 \mod p$, so we don't ever need to compute $a^2 \mod p$.
- If $a^{-1} \equiv b$, then $b^{-1} \equiv a$, so we can store two inverses each time we find an identity $ab \equiv 1 \mod p$.
- If $ab \equiv 1 \mod p$ for some $b < a$, we will have already discovered this fact when looking for $b^{-1}$. So we need only compute $ab \mod p$ for $1 < a < b < p-1$.

The brute-force algorithm suggested by the problem would perform $(p-1)^2$ multiplications modulo $p$. With the implemented optimisations, we perform at most $\frac{(p-3)^2}{4}$, and on average $\frac{(p-3)^2}{8}$.

To see that this is the case, note that in the worst case, for each $a$ we need to compute the multiplication with every $b>a$ that wasn't already known to be the inverse of some other value. So for 2, we need to compute $p-4$ multiplications; for 3, we need to compute $p-6$ multiplications, and so forth, up to $\frac{p-1}{2}$, where we perform 1 multiplication. Summing, we get $(p-3)^2/4$. For the average, we expect to find the inverse about halfway through testing values.

### Problem 2:

Write a program to turn a matrix into reduced row echelon form using Gaussian elimination. From your output, compute the ranks of each of the following matrices, and give bases for their row spaces.

$A_1 = \begin{bmatrix} 11 & 1 & 7 & 2 & 0 \\\ 8 & 0 & 2 & 5 & 11 \\\ 2 & 1 & 2 & 6 & 5 \\\ 7 & 4 & 5 & 3 & 1 \end{bmatrix} \mod 5 \text{ and} \mod 11;$  

$A_2 = \begin{bmatrix} 0 & 1 & 1 & 3 & 5 & 2 \\\ 1 & 2 & 3 & 8 & 9 & 0 \\\ 0 & 1 & 1 & 2 & 3 & 2 \\\ 2 & 1 & 3 & 7 & 9 & 1 \\\ 2 & 1 & 3 & 8 & 10 & 0 \end{bmatrix} \mod 23$

#### Solution:

This program is implemented as `helpers.gaussianElimination`. $\text{rank}(A_1)$ is 3 when working mod 5, and 4 when working mod 11; $\text{rank}(A_2) = 5$.

Considering $A_1$ mod 5, the reduced row echelon form is 

$\begin{bmatrix} 1 & 0 & 0 & 4 & 0 \\\ 0 & 1 & 0 & 0 & 4 \\\ 0 & 0 & 1 & 4 & 3 \\\ 0 & 0 & 0 & 0 & 0 \end{bmatrix}$

from which we can read off a rank of 3, and a row space basis of $\begin{bmatrix} 1 & 0 & 0 & 4 & 0 \end{bmatrix}$, $\begin{bmatrix} 0 & 1 & 0 & 0 & 4 \end{bmatrix}$, $\begin{bmatrix} 0 & 0 & 1 & 4 & 3 \end{bmatrix}$.

Considering $A_1$ mod 11, the reduced row echelon form is 

$\begin{bmatrix} 1 & 0 & 3 & 0 & 0 \\\ 0 & 1 & 7 & 0 & 0 \\\ 0 & 0 & 0 & 1 & 0 \\\ 0 & 0 & 0 & 0 & 1 \end{bmatrix}$

from which we can read off a rank of 4, and a row space basis of $\begin{bmatrix} 1 & 0 & 3 & 0 & 0 \end{bmatrix}$, $\begin{bmatrix} 0 & 1 & 7 & 0 & 0 \end{bmatrix}$, $\begin{bmatrix} 0 & 0 & 0 & 1 & 0 \end{bmatrix}$, $\begin{bmatrix} 0 & 0 & 0 & 0 & 1 \end{bmatrix}$.

Considering $A_2$ mod 23, the reduced row echelon form is 

$\begin{bmatrix} 1 & 0 & 1 & 0 & 0 & 0 \\\ 0 & 1 & 1 & 0 & 0 & 0 \\\ 0 & 0 & 0 & 1 & 0 & 0 \\\ 0 & 0 & 0 & 0 & 1 & 0 \\\ 0 & 0 & 0 & 0 & 0 & 1 \end{bmatrix}$

from which we can read off a rank of 5, and a row space basis of $\begin{bmatrix} 1 & 0 & 1 & 0 & 0 & 0 \end{bmatrix}$, $\begin{bmatrix} 0 & 1 & 1 & 0 & 0 & 0 \end{bmatrix}$, $\begin{bmatrix} 0 & 0 & 0 & 1 & 0 & 0 \end{bmatrix}$, $\begin{bmatrix} 0 & 0 & 0 & 0 & 1 & 0 \end{bmatrix}$. $\begin{bmatrix} 0 & 0 & 0 & 0 & 0 & 1 \end{bmatrix}$.

