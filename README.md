# CATAM 1.1: Matrices Over Finite Fields

A rewrite of my Cambridge CATAM undergraduate project, "Matrices Over Finite Fields".

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

### Kernel of a matrix

The [kernel](https://en.wikipedia.org/wiki/Kernel_(linear_algebra)) of a rank $r$, $m \times n$ matrix $A = (a_{ij})$ in reduced row echelon form has a simple basis. If there are pivots in columns $j_1 < j_2 < ... < j_r$, then consider $$w_k = e_k - \sum_{1 \le j_s \le k}a_{sk}e_{j_s},$$ where $e_j$ are the standard basis of $GF(p)^n$. We claim that the $w_k, k \in [n]/\{j_1, ..., j_k\}$ form a basis of the kernel. One can immediately see they are linearly independent, since the $e_j$ are, and we can compute $$(Aw_k)_i = a_{ik} - \sum_{1\le j_s\le k} a_{sk}\delta_{i, s} = a_{ik} - a_{ik} = 0,$$ where $\delta$ is the [Kronecker delta](https://en.wikipedia.org/wiki/Kronecker_delta). So the $w_k$ are a linearly independent set of size $n-r$ in the kernel of $A$, so are a basis.

When $A$ is an arbitrary matrix, that is, not in reduced row echelon form, there is some invertible matrix $R$ such that $RA$ is in reduced row echelon form. Suppose $U$ is a basis for the kernel of $RA$, derived using the above method. Then $RAU = \underline{0}$, so $AU = R^{-1}\underline{0} = \underline{0}$. So $U$ is a basis for the kernel of $A$ as well.

### Annihilators

Given a set $S\subseteq V$ of vectors, the [annihilator](https://en.wikipedia.org/wiki/Annihilator_(linear_algebra)) of $S$, written $S^{\circ}$, is the set of $f\in V^*$ such that $f(S) = \underline{0}$, where $V^*$ is the [dual space](https://en.wikipedia.org/wiki/Dual_space). If $U = \langle S\rangle$ is the [span](https://en.wikipedia.org/wiki/Linear_span) of $S$, then $U^{\circ} = S^{\circ}$, so we can always assume $S$ is a linear subspace. $S^{\circ}$ will always be a linear subspace. Moreover, $\dim U + \dim U^{\circ} = \dim V$.

In the special case that $V = F^n$, $V\cong V^*$, so both $U$ and $U^*$ can be viewed as subspaces of $V$. When $U$ is the rowspace of a matrix $A$, $U^{\circ}$ is the kernel. $U^{\circ\circ} = U$.

When $U$, $W$ are linear subspaces of $V$, we can form the subspaces $U\cap W$ and $U+W = \{u+w|u\in U, w\in W\}$. We have that $(U\cap W)^{\circ} = U^{\circ}+W^{\circ}$ and $(U+W)^{\circ} = U^{\circ}\cap W^{\circ}$.

### The isomorphism theorems for vector spaces

The [first isomorphism theorem](https://en.wikipedia.org/wiki/Isomorphism_theorems#Theorem_A_(modules)) states that, for $f:V\rightarrow W$ a linear map, $\text{Im}f \cong \frac{V}{\ker f}$, from which we can deduce that $\dim V = \dim(\text{Im}f) + \dim(\ker f)$.

The [second isomorphism theorem](https://en.wikipedia.org/wiki/Isomorphism_theorems#Theorem_B_(modules)) states that, for $U$ and $W$ subspaces of $V$, $\frac{U+W}{W} \cong \frac{U}{U\cap W}$, from which we can deduce that $\dim(U+W) = \dim U + \dim W - \dim(U\cap W)$.

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

### Problem 3:

Write a program to compute a basis for the kernel of a matrix. Describe briefly how your algorithm works. Find bases for the kernels of the matrix $A_1$ modulo 5, modulo 7 and modulo 13. Now find bases for the kernels of the matrix $A_2$ modulo every prime below 30. Do you get the same result for every prime?

#### Solution:

This program is implemented as `helpers.findKernelBasis`. The theory behind the algorithm is explained in the [Kernel of a matrix](#kernel-of-a-matrix) section.

Modulo 5, $\ker(A_1)$ has a basis of $\begin{bmatrix} 3 & 1 & 0 & 0 & 0 \end{bmatrix}^T$, $\begin{bmatrix} 2 & 0 & 1 & 1 & 0 \end{bmatrix}^T$.

Modulo 7, $\ker(A_1)$ is spanned by $\begin{bmatrix} 5 & 1 & 6 & 0 & 1 \end{bmatrix}^T$.

Modulo 13, $\ker(A_1)$ is spanned by $\begin{bmatrix} 5 & 9 & 11 & 1 & 1 \end{bmatrix}^T$.

For $A_2$, we display a table of kernel bases for each prime modulus $p<30$ below:

| $p$ | Kernel basis                                            |
| -   | ------------------------------------------------------- |
| 2   | $\begin{bmatrix} 1 & 1 & 1 & 0 & 0 & 0 \end{bmatrix}^T$ |
| 3   | $\begin{bmatrix} 2 & 2 & 1 & 0 & 0 & 0 \end{bmatrix}^T$ |
| 5   | $\begin{bmatrix} 4 & 4 & 1 & 0 & 0 & 0 \end{bmatrix}^T$ |
| 7   | $\begin{bmatrix} 6 & 6 & 1 & 0 & 0 & 0 \end{bmatrix}^T$ |
| 11  | $\begin{bmatrix} 10& 10& 1 & 0 & 0 & 0 \end{bmatrix}^T$ |
| 13  | $\begin{bmatrix} 12& 12& 1 & 0 & 0 & 0 \end{bmatrix}^T$ |
| 17  | $\begin{bmatrix} 16& 16& 1 & 0 & 0 & 0 \end{bmatrix}^T$ |
| 19  | $\begin{bmatrix} 18& 18& 1 & 0 & 0 & 0 \end{bmatrix}^T$ |
| 23  | $\begin{bmatrix} 22& 22& 1 & 0 & 0 & 0 \end{bmatrix}^T$ |
| 29  | $\begin{bmatrix} 28& 28& 1 & 0 & 0 & 0 \end{bmatrix}^T$ |

Notice that every basis is equivalent in the respective $GF(p)$ to $\begin{bmatrix} -1 & -1 & 1 & 0 & 0 & 0 \end{bmatrix}^T$. This is unsurprising: if the reduced row echelon form is the same between primes, our calculation of the basis doesn't depend on the modulus. Indeed, for any field $F$ in which $A_2$ makes sense to define and has full rank (e.g. over $\mathbb{R}$, $\mathbb{Z}$, $\mathbb{Q}$, $\mathbb{C}$, $\mathbb{Q}_p$), the kernel will have the same basis.

### Problem 4:

Use your program from Problem 3 to find $U^{\circ}$, where we work modulo 23 and $U$ is the row space of the matrix $A_1$. Similarly find $(U^{\circ})^{\circ}$ and verify that it is equal to $U$.

#### Solution:

We can compute $U^{\circ} = \ker(A_1) = \langle \begin{bmatrix} 15 & 9 & 16 & 18 & 1 \end{bmatrix}^T \rangle$, while $$U^{\circ\circ} = \langle \begin{bmatrix} 4 & 1 & 0 & 0 & 0 \end{bmatrix}, \begin{bmatrix} 2 & 0 & 1 & 0 & 0 \end{bmatrix}, \begin{bmatrix} 8 & 0 & 0 & 1 & 0 \end{bmatrix}, \begin{bmatrix} 3 & 0 & 0 & 0 & 1 \end{bmatrix} \rangle.$$ To verify that it is equal to $U$, it sufficies to express the rows of $A_1$ as sums of the basis vectors of $U^{\circ\circ}$, since $A_1$ has rank 4. Writing $e_1, e_2, e_3, e_4$ for the basis vectors in the order listed, one can check that $$A_1 = \begin{bmatrix} e_1+7e_2+2e_3 \\\ 2e_2+5e_3+11e_4 \\\ e_1+2e_2+6e_3+5e_4 \\\ 4e_1+5e_2+3e_3+e_4 \end{bmatrix} \mod 23.$$

### Problem 5:

Write a program that, given matrices $A$ and $B$ with row spaces $U$ and $W$, computes bases for $U$, $W$, $U + W$ and $U \cap W$. Explain briefly how your program works. Comment on the relationship between the dimensions of the four spaces computed. Run your program on the following examples:
- Modulo 19 with $U$ the row space of
    $$A_3 = \begin{bmatrix} 1 & 0 & 0 & 0 & 3 & 0 & 0 \\\ 
                            0 & 5 & 0 & 1 & 6 & 3 & 0 \\\
                            0 & 0 & 5 & 0 & 2 & 0 & 0 \\\
                            2 & 4 & 0 & 0 & 0 & 5 & 1 \\\
                            4 & 3 & 0 & 0 & 6 & 2 & 6 \end{bmatrix} $$
    and $W$ the kernel of $A_3$.

- Again take $U$ the row space of $A_3$ and $W$ the kernel of $A_3$, but this time modulo 7.

#### Solution:

Bases for $U$ and $W$ are given by the non-zero rows of the reduced row echelon forms of $A$ and $B$ respectively. A base for $U+W$ is derived the same way from $\begin{bmatrix} A \\\ B \end{bmatrix}$. Finally, since $U \cap W = (U\cap W)^{\circ\circ} = (U^{\circ}+W^{\circ})^{\circ}$, we can use `helpers.findKernelBasis` three times (on $A$, $B$, and $\begin{bmatrix} \ker(A)^T \\\ \ker(B)^T \end{bmatrix}$) to obtain a basis for $U \cap W$. This algorithm is implemented in `helpers.secondIsomorphismTheoremBases`.

For the given moduli and linear subspaces $U$, $W$, we obtain the following results:

| $p$ | Basis of $U$          | Basis of $W$   | Basis of $U + W$ | Basis of $U \cap W$ |
| --- | --------------------- | -------------- | ---------------- | ------------------- |
| 19  | $\begin{bmatrix} 1 & 0 & 0 & 0 & 0 \\\ 0 & 1 & 0 & 0 & 0 \\\ 0 & 0 & 1 & 0 & 0 \\\ 0 & 0 & 0 & 1 & 0 \\\ 0 & 0 & 0 & 0 & 1 \\\ 6 & 3 & 16 & 0 & 17 \\\ 1 & 14 & 9 & 8 & 6 \end{bmatrix}$ | $\begin{bmatrix} 1 & 0 \\\ 0 & 1 \\\ 9 & 0 \\\ 18 & 2 \\\ 6 & 0 \\\ 1 & 4 \\\ 12 & 14 \end{bmatrix}$ | $\begin{bmatrix} 1 & 0 & 0 & 0 & 0 & 0 & 0 \\\ 0 & 1 & 0 & 0 & 0 & 0 & 0 \\\ 0 & 0 & 1 & 0 & 0 & 0 & 0 \\\ 0 & 0 & 0 & 1 & 0 & 0 & 0 \\\ 0 & 0 & 0 & 0 & 1 & 0 & 0 \\\ 0 & 0 & 0 & 0 & 0 & 1 & 0 \\\ 0 & 0 & 0 & 0 & 0 & 0 & 1 \end{bmatrix}$ | $\begin{bmatrix} 0 \\\ 0 \\\ 0 \\\ 0 \\\ 0 \\\ 0 \\\ 0 \end{bmatrix}$ |
| 7   | $\begin{bmatrix} 1 & 0 & 0 & 0 & 0 \\\ 0 & 1 & 0 & 0 & 0 \\\ 0 & 0 & 1 & 0 & 0 \\\ 0 & 0 & 0 & 1 & 0 \\\ 0 & 0 & 0 & 0 & 1 \\\ 0 & 3 & 0 & 2 & 0 \\\ 0 & 2 & 0 & 4 & 0 \end{bmatrix}$ | $\begin{bmatrix} 0 & 0 \\\ 1 & 0 \\\ 0 & 0 \\\ 0 & 1 \\\ 0 & 0 \\\ 3 & 2 \\\ 2 & 4 \end{bmatrix}$ | $\begin{bmatrix} 1 & 0 & 0 & 0 & 0 \\\ 0 & 1 & 0 & 0 & 0 \\\ 0 & 0 & 1 & 0 & 0 \\\ 0 & 0 & 0 & 1 & 0 \\\ 0 & 0 & 0 & 0 & 1 \\\ 0 & 3 & 0 & 2 & 0 \\\ 0 & 2 & 0 & 4 & 0 \end{bmatrix}$ | $\begin{bmatrix} 0 & 0 \\\ 4 & 5 \\\ 0 & 0 \\\ 5 & 3 \\\ 0 & 0 \\\ 1 &0 \\\ 0 & 1 \end{bmatrix}$ |
