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

