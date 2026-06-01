# PFC Square Lattice Density Function: Derivation & Visualization

Derivation and visualization of the two-mode PFC density function for a square lattice, with annotated nearest-neighbor (j=1) and next-nearest-neighbor (j=2) modes.

> **Author:** Tianmeng Zhang 

---

## What this repository contains

| File | Description |
|------|-------------|
| `docs/derivation.md` | Full step-by-step mathematical derivation |
| `src/visualize_density.py` | Python visualization with annotated atom positions |
| `results/` | Generated figures (auto-created on first run) |

---

## The final result

$$n(\mathbf{r}) = A_1 \cdot \left[2\cos\frac{2\pi x}{a} + 2\cos\frac{2\pi y}{a}\right] + A_2 \cdot 4\cos\frac{2\pi x}{a}\cos\frac{2\pi y}{a}$$

---

## Derivation summary (5 steps)

**Step 1 — Fourier expansion**  
Any periodic density field can be written as $n(\mathbf{r}) = \sum_\mathbf{G} \hat{n}_\mathbf{G} e^{i\mathbf{G}\cdot\mathbf{r}}$. The complex exponential basis is chosen because each term naturally satisfies the lattice periodicity — no additional constraints needed.

**Step 2 — 8 reciprocal lattice vectors (two-mode approximation)**  
For a square lattice with $k=2\pi/a$, we keep only the two nearest shells in reciprocal space:

| Mode | Directions | $|\mathbf{G}|$ |
|------|-----------|------|
| j=1 (nearest) | $(±k,0),\ (0,±k)$ — 4 vectors | $k$ |
| j=2 (next-nearest) | $(±k,±k)$ — 4 vectors | $k\sqrt{2}$ |

**Step 3 — Pairing conjugate directions → cosines**

$$e^{+ikx} + e^{-ikx} = 2\cos(kx)$$

$$e^{+i(kx+ky)} + e^{-i(kx+ky)} + e^{+i(kx-ky)} + e^{-i(kx-ky)} = 4\cos(kx)\cos(ky)$$

The last step uses the sum-to-product identity $\cos A + \cos B = 2\cos\frac{A+B}{2}\cos\frac{A-B}{2}$.

**Step 4 — Assemble two-mode density**  
Each mode carries a free-energy-determined amplitude $A_j$, giving the final expression above.

**Step 5 — Physical verification**

| Position | j=1 | j=2 | Physical role |
|----------|-----|-----|---------------|
| Lattice site $(0,0)$ | +4 | +4 | Density maximum = **atom** |
| Edge center $(a/2,0)$ | 0 | −4 | j=2 suppresses density → **between atoms** |
| Face center $(a/2,a/2)$ | −4 | +4 | Modes partially cancel |

**j=2 sharpens atomic peaks**: it suppresses density at edge centers, making the peaks more delta-function-like and physically realistic.

---

## Figures produced

| Figure | Content |
|--------|---------|
| `fig1_reciprocal_space.png` | 8 G-vectors in reciprocal space, grouped by mode |
| `fig2_mode_j1_annotated.png` | j=1 density field with atom and edge-center annotations |
| `fig3_mode_j2_annotated.png` | j=2 density field with atom, face-center, and edge-center annotations |
| `fig4_combined_density_annotated.png` | Combined density for A₂=0, 0.2, 0.4 — showing peak sharpening |
| `fig5_1d_crosssection.png` | 1D cross-section along x-axis comparing j=1, j=2, and sum |

---

## Usage

```bash
pip install numpy matplotlib
python src/visualize_density.py
# figures saved to results/
```

---

## Connection to the PFC model

This density function is the real-space representation of the two-mode amplitude approximation used in PFC simulations. The amplitudes $A_1, A_2$ are not free parameters — they are determined by minimizing the PFC free energy functional at each temperature. When $A_1 = A_2 = 0$, the density is uniform (liquid phase); when $A_1, A_2 > 0$, the periodic density peaks appear (solid phase). The solid–liquid phase transition in PFC is precisely the point where these amplitudes grow from zero.

---

## Reference

Elder, K.R., Grant, M. (2004). Modeling elastic and plastic deformations in nonequilibrium processing using phase field crystals. *Physical Review E*, 70, 051605.
