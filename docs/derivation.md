# PFC Square Lattice Density Function: Full Derivation

**Author:** Tianmeng Zhang  
**Affiliation:** Technion – Israel Institute of Technology, MSc Materials Science & Engineering

---

## Overview

This document derives the two-mode density function used in the Phase Field Crystal (PFC) model for a square lattice, starting from the Fourier expansion of a periodic field and ending at the final expression involving nearest-neighbor (j=1) and next-nearest-neighbor (j=2) modes.

The final result is:

$$n(\mathbf{r}) = A_1 \cdot \left[2\cos\frac{2\pi x}{a} + 2\cos\frac{2\pi y}{a}\right] + A_2 \cdot 4\cos\frac{2\pi x}{a}\cos\frac{2\pi y}{a}$$

---

## Step 1: Fourier Expansion of a Periodic Density Field

Any function $n(\mathbf{r})$ satisfying the lattice periodicity condition:

$$n(\mathbf{r} + \mathbf{R}) = n(\mathbf{r}), \quad \mathbf{R} = n_1\mathbf{a}_1 + n_2\mathbf{a}_2$$

can be exactly represented as a Fourier series over reciprocal lattice vectors $\mathbf{G}$:

$$\boxed{n(\mathbf{r}) = \sum_{\mathbf{G}} \hat{n}_{\mathbf{G}}\, e^{i\mathbf{G}\cdot\mathbf{r}}}$$

**Why complex exponentials?**  
Each basis function $e^{i\mathbf{G}\cdot\mathbf{r}}$ naturally satisfies the lattice symmetry:

$$e^{i\mathbf{G}\cdot(\mathbf{r}+\mathbf{R})} = e^{i\mathbf{G}\cdot\mathbf{r}} \cdot \underbrace{e^{i\mathbf{G}\cdot\mathbf{R}}}_{=1} = e^{i\mathbf{G}\cdot\mathbf{r}}$$

The condition $e^{i\mathbf{G}\cdot\mathbf{R}} = 1$ is exactly the definition of reciprocal lattice vectors. Every term in the expansion is individually periodic — no additional constraints needed.

**The density field $n(\mathbf{r})$ is real**, which requires $\hat{n}_{-\mathbf{G}} = \hat{n}^*_{\mathbf{G}}$. This means positive and negative $\mathbf{G}$ directions always appear in conjugate pairs, and their imaginary parts cancel:

$$\hat{n}_{\mathbf{G}} e^{i\mathbf{G}\cdot\mathbf{r}} + \hat{n}_{-\mathbf{G}} e^{-i\mathbf{G}\cdot\mathbf{r}} = 2\,\text{Re}[\hat{n}_{\mathbf{G}} e^{i\mathbf{G}\cdot\mathbf{r}}]$$

The complex exponential is a computational convenience, not a physical imaginary quantity.

---

## Step 2: Reciprocal Lattice Vectors of the Square Lattice

For a square lattice with lattice constant $a$, define $k = 2\pi/a$.

Reciprocal lattice vectors take the form $\mathbf{G} = k(h, l)$ for integers $h, l$.  
We group all $\mathbf{G}$ vectors by their magnitude $|\mathbf{G}|$:

### Mode j=1 — Nearest neighbors in reciprocal space

$$|\mathbf{G}| = k\sqrt{1^2+0^2} = k = \frac{2\pi}{a}$$

| Direction index $l$ | $\mathbf{G}_l$ | $\mathbf{G}_l \cdot \mathbf{r}$ | Exponential term |
|---|---|---|---|
| 1 | $(+k,\; 0)$ | $+kx$ | $e^{+ikx}$ |
| 2 | $(-k,\; 0)$ | $-kx$ | $e^{-ikx}$ |
| 3 | $(0,\; +k)$ | $+ky$ | $e^{+iky}$ |
| 4 | $(0,\; -k)$ | $-ky$ | $e^{-iky}$ |

### Mode j=2 — Next-nearest neighbors in reciprocal space

$$|\mathbf{G}| = k\sqrt{1^2+1^2} = k\sqrt{2} = \frac{2\pi\sqrt{2}}{a}$$

| Direction index $l$ | $\mathbf{G}_l$ | $\mathbf{G}_l \cdot \mathbf{r}$ | Exponential term |
|---|---|---|---|
| 5 | $(+k,\; +k)$ | $+kx+ky$ | $e^{+i(kx+ky)}$ |
| 6 | $(+k,\; -k)$ | $+kx-ky$ | $e^{+i(kx-ky)}$ |
| 7 | $(-k,\; +k)$ | $-kx+ky$ | $e^{-i(kx-ky)}$ |
| 8 | $(-k,\; -k)$ | $-kx-ky$ | $e^{-i(kx+ky)}$ |

**Why these two modes?**  
The 4-fold rotational symmetry of the square lattice makes all directions within the same mode physically equivalent — they all carry the same amplitude $A_j$. Modes beyond j=2 contribute negligibly to the free energy (their amplitudes decay rapidly), so the two-mode approximation captures the essential physics.

---

## Step 3: Pairing Conjugate Directions — Eliminating Imaginary Parts

### j=1 mode (4 directions → 2 cosines)

Apply Euler's formula $e^{i\theta} + e^{-i\theta} = 2\cos\theta$ to each conjugate pair:

**Pair (1,2) — x-direction:**

$$e^{+ikx} + e^{-ikx} = 2\cos(kx) = 2\cos\!\frac{2\pi x}{a}$$

**Pair (3,4) — y-direction:**

$$e^{+iky} + e^{-iky} = 2\cos(ky) = 2\cos\!\frac{2\pi y}{a}$$

Sum of all j=1 terms:

$$\boxed{\sum_{l=1}^{4} e^{i\mathbf{G}_l \cdot \mathbf{r}} = 2\cos\!\frac{2\pi x}{a} + 2\cos\!\frac{2\pi y}{a}}$$

### j=2 mode (4 directions → 1 product of cosines)

**Pair (5,8) — diagonal (+,+) and (−,−):**

$$e^{+i(kx+ky)} + e^{-i(kx+ky)} = 2\cos(kx+ky)$$

**Pair (6,7) — diagonal (+,−) and (−,+):**

$$e^{+i(kx-ky)} + e^{-i(kx-ky)} = 2\cos(kx-ky)$$

**Sum of all four j=2 terms:**

$$2\cos(kx+ky) + 2\cos(kx-ky)$$

Apply the sum-to-product identity $\cos A + \cos B = 2\cos\!\frac{A+B}{2}\cos\!\frac{A-B}{2}$:

$$= 2 \cdot 2\cos\!\left(\frac{(kx+ky)+(kx-ky)}{2}\right)\cos\!\left(\frac{(kx+ky)-(kx-ky)}{2}\right)$$

$$= 4\cos(kx)\cos(ky)$$

$$\boxed{\sum_{l=5}^{8} e^{i\mathbf{G}_l \cdot \mathbf{r}} = 4\cos\!\frac{2\pi x}{a}\cos\!\frac{2\pi y}{a}}$$

---

## Step 4: Assemble the Two-Mode Density Function

Each mode $j$ carries an amplitude $A_j$ (determined by free energy minimization). Combining both modes:

$$\boxed{n(\mathbf{r}) = A_1\!\left[2\cos\!\frac{2\pi x}{a} + 2\cos\!\frac{2\pi y}{a}\right] + A_2 \cdot 4\cos\!\frac{2\pi x}{a}\cos\!\frac{2\pi y}{a}}$$

---

## Step 5: Physical Verification at Special Points

Set $a=1$ for simplicity. Evaluate $n$ at three characteristic positions:

| Position | $n_{j=1}$ | $n_{j=2}$ | Role |
|---|---|---|---|
| Lattice site $(0,0)$ | $2(+1)+2(+1)=+4$ | $4(+1)(+1)=+4$ | Both modes add → **density maximum = atom** |
| Edge center $(½,0)$ | $2(-1)+2(+1)=0$ | $4(-1)(+1)=-4$ | j=2 pushes density negative → **between atoms** |
| Face center $(½,½)$ | $2(-1)+2(-1)=-4$ | $4(-1)(-1)=+4$ | Modes partially cancel |

**The role of j=2 is to sharpen the atomic peaks**: j=1 alone leaves residual positive density at edge centers; j=2 subtracts it, making the atomic density peak more localized and delta-function-like.

---

## Physical Meaning Summary

| Symbol | Meaning |
|---|---|
| $n(\mathbf{r})$ | Time-averaged single-particle number density at position $\mathbf{r}$ |
| $a$ | Lattice constant of the square phase |
| $k = 2\pi/a$ | Magnitude of first reciprocal lattice vector |
| $A_1$ | Amplitude of j=1 mode; $A_1=0$ → liquid, $A_1>0$ → solid |
| $A_2$ | Amplitude of j=2 mode; sharpens atomic peaks |
| Density maximum | Location of atoms (lattice sites) |
| $A_j \to 0$ | Solid–liquid transition in PFC model |

---

## Reference

Elder, K.R., Grant, M. (2004). Modeling elastic and plastic deformations in nonequilibrium processing using phase field crystals. *Physical Review E*, 70, 051605.
