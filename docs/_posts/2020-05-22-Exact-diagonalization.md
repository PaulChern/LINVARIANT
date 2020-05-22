---
layout: post
title:  "Exact diagonalization: Bose–Hubbard model"
date:   2020-05-22 01:30:13 +0800
categories: Physics
tags: tutorial
comments: 1
lang: en-US
math: true
description: To calculate the ground state and the corresponding energy of a quantum system is a very important task in physics. As we all known, the ground state energy is just the smallest eigenvalue of the Hamiltonian of this system.
---

## Introduction

To calculate the **ground state** and the corresponding energy of a quantum system is a very important task in physics. As we all known, the ground state energy is just the smallest eigenvalue of the Hamiltonian of this system. Therefore, we can directly find the ground state and its energy once given the Hamiltonian in principle.

$$
\gdef\ee{\mathrm{e}}
\gdef\ii{\mathrm{i}}
\gdef\bra#1{\langle#1\vert}
\gdef\ket#1{\vert#1\rangle}
\gdef\ev#1{\langle#1\rangle}
$$

For a specific set of basis $\ket{\psi_i}$ (a basis transformation), the Hamiltonian operator $\hat{H}$ can be expressed with a matrix:

$$
H_{ij} = \bra{\psi_i}\hat{H}\ket{\psi_j}.
$$

From the knowledge of linear algebra, we know that the eigenvalues of a matrix will not change under such a basis transformation (namely similarity transformation). So the ground state energy is simply the smallest eigenvalue of the matrix $\bm{H}$. To get the eigenvalue, there are already some powerful numerical algorithms, such as the [**Lanczos algorithm**](https://en.wikipedia.org/wiki/Lanczos_algorithm).

Now the problem becomes how to find a "good" set of basis such that $\bm{H}$ can have a simple form. If one is lucky enough that $\bm{H}$ is diagonal under his basis, then there is no need to do more thing since he has already got the solution. Clearly, it's impossible to do that in reality, but with the help of symmetry of the system, we can still greatly simplify the matrix (into block diagonal of band matrix, for example).

## Bose--Hubbard model

Now let's turn to a concrete example. The Hamiltonian of the one dimensional [**Bose--Hubbard model**](https://en.wikipedia.org/wiki/Bose%E2%80%93Hubbard_model) can be written as:

$$
\hat{H} = -t \sum_{\ev{ij}} \left(\hat{a}_i^\dagger\hat{a}_j + \hat{a}_j^\dagger\hat{a}_i\right)
          + \frac{U}{2} \sum_{i=1}^N  \,\hat{n}_i \left(\hat{n}_i - 1\right)
          - \mu \sum_{i=1}^N  \,\hat{n}_i,
$$

where $\hat{a}_i^\dagger$ and $\hat{a}_i$ are creation and annihilation operator respectively (recall the harmonical oscillator), $\hat{n}_i=\hat{a}_i^\dagger\hat{a}_i$ is the particle number operator, $N$ is the total numbers of the sites (not the particles!).

In the first part, $\hat{a}_i^\dagger\hat{a}_j$ means that one boson is annihilated at site $j$ and then created at site $i$, so the net effect is this boson "moves" from $j$ to $i$. That's why we call it **hopping term**. Summation indices $\ev{ij}$ means that only neighboring lattice sites are considered, since the particles can't move too far.

The second part is the **interaction term**, which describes the on-site interaction of the boson. It has the same form as the electrostatic energy of $N$ electric charges. In the following article, we assume $U>0$ so that the interaction is repulsive.

The **chemical potential** $\mu$ in the last part essentially sets the total number of particles. If the number is fixed, then it will become a constant and we can simply ignore it.

To represent the wave fuction, we can choose the occupation number basis:

$$
\ket{n_1,\, n_2, \, \ldots,\, n_i,\, \ldots,\, n_N}.
$$

Here, $N$ is the number of sites, while $M=\sum_{i=1}^N n_i$ is the total number of particles. Apply the operators on it, we have:

$$
\begin{aligned}
   \vphantom{\hat{a}_i^\dagger}
   \hat{a}_i         \, \ket{\ldots,\, n_i  ,\, \ldots}
&= \sqrt{n_i}        \, \ket{\ldots,\, n_i-1,\, \ldots}, \\
   \hat{a}_i^\dagger \, \ket{\ldots,\, n_i  ,\, \ldots}
&= \sqrt{n_i+1}      \, \ket{\ldots,\, n_i+1,\, \ldots}, \\
   \vphantom{\hat{a}_i^\dagger}
   \hat{n}_i         \, \ket{\ldots,\, n_i  ,\, \ldots}
&= n_i               \, \ket{\ldots,\, n_i  ,\, \ldots}.
\end{aligned}
$$

As expected, the occupation number basis is the eigenstate of $\hat{n}_i$, with eigenvalue $n_i$.

## Almost trivial solutions

When the lattice is small, we can find the ground state even by hand. Now suppose there are only two sites and two particles. All the possible configurations are the following:

$$
\ket{2,\,0}, \quad \ket{1,\,1} \quad \text{and} \quad \ket{0,\,2}.
$$

The Hamiltonian becomes

$$
\hat{H} = -t \left(\hat{a}_1^\dagger\hat{a}_2 + \hat{a}_2^\dagger\hat{a}_1\right)
          + \frac{U}{2} \Bigl[  \hat{n}_1 \left(\hat{n}_1-1\right)
                              + \hat{n}_2 \left(\hat{n}_2-1\right) \Bigr],
$$

then

$$
\begin{aligned}
\hat{H} \, \ket{2,\,0} &= -t \left(0 + \sqrt{2} \, \hat{a}_2^\dagger \, \ket{1,\,0} \right)
                        + \frac{U}{2} \Bigl[2\times(2-1) + 0\Bigr] \, \ket{2,\,0} \\
                       &= U \, \ket{2,\,0} - \sqrt{2} \, t \, \ket{1,\,1},
\end{aligned}
$$

and similarly

$$
\begin{aligned}
\hat{H} \, \ket{1,\,1} &= -\sqrt{2} \, t \, \ket{2,\,0} - \sqrt{2} \, t \, \ket{0,\,2}, \\
\hat{H} \, \ket{0,\,2} &= -\sqrt{2} \, t \, \ket{1,\,1} + U \, \ket{0,\,2}.
\end{aligned}
$$

Let's rewrite the above basis as vectors (we can see here that the Hilbert space is 3-dimensional):

$$
\begin{bmatrix} 1 \\ 0 \\ 0 \end{bmatrix} \coloneqq \ket{2,\,0}, \quad
\begin{bmatrix} 0 \\ 1 \\ 0 \end{bmatrix} \coloneqq \ket{1,\,1}, \quad
\begin{bmatrix} 0 \\ 0 \\ 1 \end{bmatrix} \coloneqq \ket{0,\,2},
$$

then the Hamiltonian can be written as a matrix:

$$
\bm{H} =
  \begin{bmatrix}
    U            & -\sqrt{2}\,t & 0            \\
    -\sqrt{2}\,t & 0            & -\sqrt{2}\,t \\
    0            & -\sqrt{2}\,t & U
  \end{bmatrix}.
$$

We can evaluate the eigensystem directly by Mathematica:

```_wl
{% raw %}H = {{U, -Sqrt[2] t, 0}, {-Sqrt[2] t, 0, -Sqrt[2] t}, {0, -Sqrt[2] t, U}};
FullSimplify @ Eigensystem @ H{% endraw %}
```

The result is

$$
\begin{aligned}
E_1    &= \frac{U - \sqrt{16t^2 + U^2}}{2}, &
\psi_1 &= \begin{bmatrix} 2\sqrt{2}\,t \\ U + \sqrt{16t^2+U^2} \\ 2\sqrt{2}\,t \end{bmatrix}; \\
E_2    &= U, &
\psi_2 &= \begin{bmatrix} 1 \\ 0 \\ -1 \end{bmatrix}; \\
E_3    &= \frac{U + \sqrt{16t^2 + U^2}}{2}, &
\psi_3 &= \begin{bmatrix} 2\sqrt{2}\,t \\ U - \sqrt{16t^2+U^2} \\ 2\sqrt{2}\,t \end{bmatrix}.
\end{aligned}
$$

Clearly, $E_1$ is the ground state energy and the state can be written as

$$
\ket{\psi_\text{ground}} = C \left[
  2\sqrt{2}\,t \Bigl(\ket{2,\,0} + \ket{0,\,2}\Bigr) +
  \left(U + \sqrt{16t^2+U^2}\right) \ket{1,\,1}
\right],
$$

where $C$ is a normalization constant.

## Symmetries of the Hamiltonian

The Hamiltonian of Bose--Hubbard model possesses a **global *U*(1) symmetry**. With

$$
\hat{a}_i \to \ee^{\ii\theta}\hat{a}_i, \quad
\hat{a}_i^\dagger \to \ee^{-\ii\theta}\hat{a}_i^\dagger \quad \text{and hence} \quad
\hat{n}_i \to \hat{n}_i,
$$

we have (ignore the chemical potential)

$$
\begin{aligned}
\hat{H} &= -t \sum_{\ev{ij}}
              \left(\hat{a}_i^\dagger\hat{a}_j + \hat{a}_j^\dagger\hat{a}_i\right)
           + \frac{U}{2} \sum_{i=1}^N  \,\hat{n}_i \left(\hat{n}_i - 1\right) \\
        & \to -t \sum_{\ev{ij}}
              \left(\ee^{-\ii\theta}\hat{a}_i^\dagger \cdot \ee^{\ii\theta}\hat{a}_j +
                    \ee^{-\ii\theta}\hat{a}_j^\dagger \cdot \ee^{\ii\theta}\hat{a}_i\right)
           + \frac{U}{2} \sum_{i=1}^N  \,\hat{n}_i \left(\hat{n}_i - 1\right)
        = \hat{H}.
\end{aligned}
$$

If we use the periodic boundary condition, there will be a **translation symmetry**. Take $i \to i+1$, the summations will be invariant:

$$
\begin{aligned}
\sum_{i=1}^N f[i] \to \sum_{i=1}^N f[i+1] = \sum_{i=0}^{N-1} f[i]
  &= f[0] + \sum_{i=1}^{N-1} f[i] \\
  &= f[N] + \sum_{i=1}^{N-1} f[i] = \sum_{i=1}^N f[i], \\
\sum_{\ev{ij}}^N g[i,\,j] \to \sum_{\ev{ij}}^N g[i+1,\,j+1]
  &= \sum_{i=1}^N \underbrace{\sum_{|j-i|=1} g[i+1,\,j+1]}_{f[i+1]} \\
  &= \sum_{i=1}^N \sum_{|j-i|=1} g[i,\,j] = \sum_{\ev{ij}}^N g[i,\,j],
\end{aligned}
$$

and so is the Hamiltonian. We can imagine the one-dimensional boson chain as a circle, with the $N$ sites placed equidistantly. Obviously the system will be invariant if rotate the circle by one unit.

This circle, or more precisely, a regular $N$-sided polygon, has another **reflection symmetry**, which means it's invariant under $i \to N-i$. This can be easily proved with the same procedure as above.

The translation and reflection symmetry make up a [**dihedral group**](https://en.wikipedia.org/wiki/Dihedral_group) $D_N$. Therefore, the total symmetry group of the one dimensional Bose--Hubbard model is $U(1) \otimes D_N$.

## Simulation

The general algorithm, as we have introduced above, is as following:

1. Choose a set basis.
2. Write down the Hamiltonian matrix in this basis.
3. Use Lanczos algorithm to find the eigenvalues.

Now let's use Mathematica to turn this idea into real-world code.

### Basis

A basis vector is a `List` containing `siteNum` elements, whose sum is `particleNum`. This job can be done via the build-in function `IntegerPartitions`. Then we use `Permutations` to get all the possible configurations of these particles. Finally we use `Flatten` and `ReverseSort` to make the result prettier:

```_wl
getBasis[siteNum_, particleNum_] :=
  ReverseSort @ Catenate[
    Permutations[PadRight[#, siteNum]] & /@ IntegerPartitions[particleNum, siteNum]]
```

Here is a simple example (you can easily check it by hand):

```_wl
In[1] := getBasis[3, 3] // TableForm
Out[1]//TableForm =
    3    0    0
    2    1    0
    2    0    1
    1    2    0
    1    1    1
    1    0    2
    0    3    0
    0    2    1
    0    1    2
    0    0    3
```

A more reliable test is to calculate the dimension of the Hilbert space when `siteNum` equals `particleNum`:

```_wl
In[2] := Length[getBasis[#, #]] & /@ Range[10]
Out[2] = {1, 3, 10, 35, 126, 462, 1716, 6435, 24310, 92378}
```

It is actually the [A001700](https://oeis.org/A001700) sequence. In general case, the dimension is

$$
D_{N,\,M} = \text{length of \texttt{getBasis[N, M]}}
          = \binom{N+M-1}{N} = \frac{(N+M-1)!}{N! \, (M-1)!}.
$$

### Matrix

First, let's write down the Hamiltonian again for later reference:

$$
\hat{H} = - \sum_{\ev{ij}} \left(\hat{a}_i^\dagger\hat{a}_j + \hat{a}_j^\dagger\hat{a}_i\right)
          + \frac{U}{2t} \sum_{i=1}^N  \,\hat{n}_i \left(\hat{n}_i - 1\right).
$$

Note that the constants have been merged for simplicity.

As denoted before, the occupation number basis is the eigenstate of $\hat{n}_i$, so the second part is already diagonalized. The first part is off-diagonal. But fortunately, it consists merely adjacent elements, hence the matrix is still sparse. We formally write the matrix as a combination of two parts:

```_wl
getMatrix[basis_, couplingConst_] :=
  With[{basisNumRange = Range @ Length @ basis},
    SparseArray @ Join[
      kineticPart[basis, AssociationThread[basis -> basisNumRange], basisNumRange],
      interactionPart[basis, couplingConst, basisNumRange]]]
```

Here, `basis` is a list of vectors generated by `getBasis[]` and `couplingConst` equals to $U/t$. The `interactionPart[]` is straightforward:

```_wl
interactionPart[basis_, couplingConst_, basisNumRange_] :=
  MapThread[{#1, #1} -> #2 &,
    {basisNumRange, 0.5 * couplingConst * Sum[i * (i - 1), {i, #}] & /@ basis}]
```

`interactionPart[]` returns a list of rules like `{i, j} -> value`. They are taken by `SparseArray[]` to construct the matrix then.

The `kineticPart[]` (i.e. hopping term), however, is much more complicated and need some explanation:

```_wl
{% raw %}kineticPart[basis_, positionMap_, basisNumRange_] :=
  Catenate @ MapThread[kineticPartMapFunc] @ {
    Apply[{positionMap[#1], #2} &, DeleteCases[{_, 0.}] /@
      Transpose[{opADagAState[basis], opADagAValue[basis]}, {3, 1, 2}], {2}],
    basisNumRange}
opADagAState[basis_] :=
  With[{len = Length @ First @ basis},
    Outer[Plus, basis, #, 1] & @ Catenate[
      NestList[RotateRight, PadRight[#, len], len - 1] & /@ {{1, -1}, {-1, 1}}]]
opADagAValue[basis_] :=
  -Sqrt[(#1 + 1.) * #2] & @@@ (Join[#, Reverse[#, {2}]] & @ Partition[#, 2, 1, 1]) & /@
    basis
kineticPartMapFunc[stateValuePairs_, index_] :=
  ({index, #1} -> #2) & @@@ stateValuePairs{% endraw %}
```

The basic idea is to calculate all the non-vanishing matrix elements. I use some list tricks to find the components of $\hat{a}_i^\dagger\hat{a}_j$ in `opADagAState[]`, then obtain their indices from the association `positionMap`. In `opADagAValue[]`, the coefficients are calculated. Note that the zero elements will be removed in `kineticPart[]`, as they make no contributions to the Hamiltonian matrix.

<figure>
  <img src="/images/exact-diagonalization/hamiltonian-matrix-n=3.svg" alt="hamiltonian-matrix-n=3" style="width: 25%;" class="invert"><img src="/images/exact-diagonalization/hamiltonian-matrix-n=4.svg" alt="hamiltonian-matrix-n=4" style="width: 25%;" class="invert"><img src="/images/exact-diagonalization/hamiltonian-matrix-n=5.svg" alt="hamiltonian-matrix-n=5" style="width: 25%;" class="invert"><img src="/images/exact-diagonalization/hamiltonian-matrix-n=6.svg" alt="hamiltonian-matrix-n=6" style="width: 25%;" class="invert">
  <figcaption markdown="span">The result of `getMatrix[]` with `siteNum = particleNum = 3, 4, 5, 6` and `couplingConst = 1`. Colors are not scaled equally in above 4 figures.</figcaption>
</figure>

### Eigensystem

The Lanczos algorithm is called `Arnoldi` in Mathematica. The implementation is based on the [ARPACK library](https://en.wikipedia.org/wiki/ARPACK) and is most useful for large sparse matrices[^wolfram-eigenvalues]. Note that we need the smallest eigenvalue for the ground state, but `Eigensystem[]` can only find the largest ones. So there are two extra minus signs here:

[^wolfram-eigenvalues]: [Eigenvalues&mdash;Wolfram Language Documentation](https://reference.wolfram.com/language/ref/Eigenvalues.html)

```_wl
getEigensystem[matrix_, num_: 1] :=
  -Eigensystem[N[-matrix], num, Method -> {"Arnoldi", "Criteria" -> "RealPart"}]
getGroundEigensystem[matrix_] := First /@ getEigensystem[matrix, 1]
```

## Numeric results

Once the ground state is obtained, we can then calculate some physical quantities, such as:

- Single-particle density matrix (SPDM)
- Condensate fraction
- Occupation variance

### SPDM

The single-particle density matrix (SPDM) is defined as

$$
\rho_{ij}^{(1)} = \bra{\psi_0} \hat{a}_i^\dagger\hat{a}_j \ket{\psi_0},
$$

where $\ket{\psi_0}$ is the ground state.

Here, we first calculate the result of $\hat{a}_i\ket{\psi_0}$, then merge and sum over the non-vanishing terms.
We only need to consider the off-diagonal elements, since $\rho\_{ii}^{(1)}$ are automatically zero.

```_wl
{% raw %}getSPDM[groundState_, basis_, {i_, i_}] := 0.
getSPDM[groundState_, basis_, {i_, j_}] :=
  Total @ Merge[AssociationThread /@ {
      basis -> groundState,
      (basis + ConstantArray[ReplacePart[
          ConstantArray[0, Length @ First @ basis], {i -> 1, j -> -1}],
          Length @ basis])
        -> (Sqrt[(#[[i]] + 1.) * #[[j]]] & /@ basis) * groundState
    }, getSPDMMergeFunc]
getSPDMMergeFunc[{a_}] := 0.
getSPDMMergeFunc[{a_, b_}] := a * b{% endraw %}
```

### Condensate fraction

The condensate fraction is essentially the largest eigenvalue $\lambda_1$ of SPDM $\bm{\rho}^{(1)}$:

$$
f_{\mathrm{c}} = \frac{\lambda_1}{M}.
$$

Note that $\rho_{ij}^{(1)}$ has the following symmetries:

$$
\rho_{ij}^{(1)} = \rho_{ji}^{(1)}, \quad \rho_{ij}^{(1)} = \rho_{i+1,\,j+1}^{(1)},
$$

as well as the periodic boundary condition. So we just need to evaluate about $\lceil(N+1)/2\rceil$ elements:

```_wl
{% raw %}getCondensateFraction[groundState_, basis_, siteNum_, particleNum_] :=
  With[{eval = Table[getSPDM[groundState, basis, {1, i}],
      {i, 2, Ceiling[(siteNum + 1) / 2]}]},
    First[Eigenvalues[#, 1]] / particleNum & @
      NestList[RotateRight,
        Join[{1.}, eval, Reverse[eval][[2 - Mod[siteNum, 2] ;;]]], siteNum - 1]]{% endraw %}
```

As its name indicated, the condensate fraction is a signal of phase transition. If $f_{\mathrm{c}} \sim 1$, we say that the system is in a condensate state. This condensate is also associated with the **off-diagonal long range order**, which can be shown with the off-diagonal elements in $\rho_{ij}^{(1)}$.

### Occupation variance

The occupation variance, or the fluctuation of the occupation number, is

$$
\sigma_i = \sqrt{\bra{\psi_0} \hat{n}_i^2 \ket{\psi_0} - \bra{\psi_0} \hat{n}_i \ket{\psi_0}^2}
$$

It can be calculated without too many tricks:

```_wl
{% raw %}getOccupationVariance[groundState_, basis_, i_: 1] := With[{
    groundStateSq = groundState * groundState, basisAtI = basis[[All, i]]},
  Sqrt[groundStateSq . (basisAtI * basisAtI) - (groundStateSq . basisAtI)^2]]{% endraw %}
```

### Plots

The above numeric results are the following:

<figure>
  <img src="/images/exact-diagonalization/numeric-results.svg" alt="numeric-results" class="invert">
  <figcaption>The SPDM is actually $\rho_{1,\,\lceil n/2\rceil}^{(1)}$. The $x$-axis is coupling constant $U/t$. Different curves are for different $M$ and $N$ (we take $M=N$ here for convenience).</figcaption>
</figure>

We can find the condensate at around $U/t=8$ clearly. SPDM and condensate fraction suffer from a strong finite size effect, while the occupation variance is size-insensitive.

In terms of performance, if we choose $M=N=12$, it takes about 240\,s to set up the matrix, 80\,s to find the ground state, 20\,s for SPDM and 90\,s for condensate fraction. It can be ignored (less than 1\,s) for other operations. The time consuming increases rapidly as the lattice become larger, so even for $N=20$ it will take too much time hence can't be done on a PC. The above test is on a MacBook Pro with a 2.3\,GHz Intel Core i5 processor.

## Note

This article is based on arXiv:1102.4006[^arxiv-1102-4006]. Here we mainly focused on the implementation of the algorithm. More details, especially physical explanation and derivations, can be found in that paper.

[^arxiv-1102-4006]: J M Zhang and R X Dong. *Exact diagonalization: the Bose--Hubbard model as an example*, [arXiv:1102.4006](https://arxiv.org/abs/1102.4006)

The complete Mathematica code can be found in [stone-zeng/toys](https://github.com/stone-zeng/toys/blob/master/exact-diagonalization/bose-hubbard.wl).

## References

<div id="footnotes"></div>
