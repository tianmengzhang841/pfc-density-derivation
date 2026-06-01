"""
PFC Square Lattice Density Function: Visualization
===================================================
Visualizes the two-mode PFC density function for a square lattice,
with explicit annotation of nearest-neighbor (j=1) and
next-nearest-neighbor (j=2) atomic positions.

Produces five figures:
  1. Real space lattice: NN and NNN labeled with arrows and wave functions
  2. j=1 mode density field with atom annotations
  3. j=2 mode density field with atom annotations
  4. Combined density n(r) = A1*j1 + A2*j2 for three values of A2
  5. 1D cross-section along x-axis

Author: Tianmeng Zhang
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch, FancyArrow
import os

# ── Parameters ────────────────────────────────────────────────────────────────
a  = 1.0
k  = 2 * np.pi / a
A1 = 0.5
A2 = 0.2
N  = 400
SAVE_DIR = "results"
os.makedirs(SAVE_DIR, exist_ok=True)

# ── Grid ──────────────────────────────────────────────────────────────────────
x = np.linspace(0, 2*a, N)
y = np.linspace(0, 2*a, N)
X, Y = np.meshgrid(x, y)

def mode_j1(X, Y, k):
    return 2*np.cos(k*X) + 2*np.cos(k*Y)

def mode_j2(X, Y, k):
    return 4*np.cos(k*X) * np.cos(k*Y)

def density(X, Y, k, A1, A2):
    return A1 * mode_j1(X, Y, k) + A2 * mode_j2(X, Y, k)

J1      = mode_j1(X, Y, k)
J2      = mode_j2(X, Y, k)
N_field = density(X, Y, k, A1, A2)

def lattice_sites(a, n_cells=2):
    sites = []
    for ix in range(n_cells + 1):
        for iy in range(n_cells + 1):
            sites.append((ix * a, iy * a))
    return sites

def face_centers(a, n_cells=2):
    sites = []
    for ix in range(n_cells):
        for iy in range(n_cells):
            sites.append((ix*a + a/2, iy*a + a/2))
    return sites

def edge_centers(a, n_cells=2):
    sites = []
    for ix in range(n_cells):
        for iy in range(n_cells + 1):
            sites.append((ix*a + a/2, iy*a))
            sites.append((iy*a, ix*a + a/2))
    return list(set(sites))

ATOM_SITES = lattice_sites(a, n_cells=2)
FACE_SITES = face_centers(a, n_cells=2)
EDGE_SITES = [(s[0], s[1]) for s in edge_centers(a, n_cells=2)
              if 0 <= s[0] <= 2*a and 0 <= s[1] <= 2*a]

def add_colorbar(fig, ax, im, label):
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label(label, fontsize=10)
    return cb

# ═══════════════════════════════════════════════════════════════════════════════
# Figure 1: Real space lattice — NN and NNN labeled
# ═══════════════════════════════════════════════════════════════════════════════
fig1, ax1 = plt.subplots(figsize=(8, 7))
ax1.set_aspect('equal')
ax1.set_xlim(-1.8*a, 1.8*a)
ax1.set_ylim(-1.8*a, 2.2*a)
ax1.axis('off')

# background grid atoms
bg_positions = [(-a,-a),(-a,0),(-a,a),(0,-a),(0,a),(a,-a),(a,a)]
for (bx, by) in bg_positions:
    ax1.add_patch(plt.Circle((bx,by), 0.10, color='lightgray',
                              ec='gray', lw=0.8, zorder=1))

# NNN atoms (amber) — 4 diagonal neighbors
nnn_positions = [(-a,-a),( a,-a),(-a, a),( a, a)]
nnn_labels    = ['(−a,−a)', '(+a,−a)', '(−a,+a)', '(+a,+a)']
nnn_wf        = ['e^{−i(kx+ky)}', 'e^{+i(kx−ky)}', 'e^{−i(kx−ky)}', 'e^{+i(kx+ky)}']
for (nx, ny) in nnn_positions:
    ax1.add_patch(plt.Circle((nx, ny), 0.18, color='#FAEEDA',
                              ec='#854F0B', lw=2, zorder=3))
    ax1.text(nx, ny, 'NNN', ha='center', va='center',
             fontsize=7, color='#633806', fontweight='bold', zorder=4)

# NN atoms (blue) — 4 axis-aligned neighbors
nn_positions = [(0, a), (0,-a), (-a, 0), (a, 0)]
nn_labels    = ['(0,+a)', '(0,−a)', '(−a,0)', '(+a,0)']
nn_wf        = ['cos(ky)', 'cos(ky)', 'cos(kx)', 'cos(kx)']
for (nx, ny) in nn_positions:
    ax1.add_patch(plt.Circle((nx, ny), 0.18, color='#E6F1FB',
                              ec='#185FA5', lw=2, zorder=3))
    ax1.text(nx, ny, 'NN', ha='center', va='center',
             fontsize=8, color='#0C447C', fontweight='bold', zorder=4)

# central atom (purple)
ax1.add_patch(plt.Circle((0, 0), 0.22, color='#534AB7',
                          ec='#3C3489', lw=2.5, zorder=5))
ax1.text(0, 0, '0', ha='center', va='center',
         fontsize=11, color='white', fontweight='bold', zorder=6)

# NN solid arrows
arrowkw_nn = dict(color='#185FA5', lw=2, length_includes_head=True,
                  head_width=0.05, head_length=0.07)
for (nx, ny) in nn_positions:
    dx, dy = nx*0.75, ny*0.75
    ax1.arrow(0, 0, dx, dy, **arrowkw_nn, zorder=2)

# NNN dashed arrows
for (nx, ny) in nnn_positions:
    dx, dy = nx*0.75, ny*0.75
    ax1.annotate('', xy=(dx, dy), xytext=(0, 0),
                 arrowprops=dict(arrowstyle='->', color='#BA7517',
                                 lw=1.8, linestyle='dashed'),
                 zorder=2)

# distance annotations
ax1.annotate('', xy=(a, 0.05), xytext=(0, 0.05),
             arrowprops=dict(arrowstyle='<->', color='#185FA5', lw=1.2))
ax1.text(a/2, 0.12, 'a', ha='center', color='#185FA5', fontsize=11)

ax1.annotate('', xy=(a*0.75, a*0.75+0.06), xytext=(0.06, 0.06),
             arrowprops=dict(arrowstyle='<->', color='#BA7517',
                             lw=1.2, linestyle='dashed'))
ax1.text(a*0.46, a*0.54, r'$a\sqrt{2}$', color='#BA7517', fontsize=10)

# wave function labels next to each NN/NNN
for (nx, ny), lbl, wf in zip(nn_positions, nn_labels, nn_wf):
    off_x = 0.25 if nx >= 0 else -0.25
    off_y = 0.25 if ny >= 0 else -0.25
    if nx == 0: off_x = 0.28
    if ny == 0: off_y = 0.22
    ax1.text(nx + off_x, ny + off_y,
             f'{lbl}\n→ {wf}',
             ha='center', va='center', fontsize=8,
             color='#0C447C',
             bbox=dict(boxstyle='round,pad=0.2', fc='#E6F1FB',
                       ec='#185FA5', alpha=0.85))

for (nx, ny), lbl, wf in zip(nnn_positions, nnn_labels, nnn_wf):
    off_x = 0.32 if nx > 0 else -0.32
    off_y = 0.32 if ny > 0 else -0.32
    ax1.text(nx + off_x, ny + off_y,
             f'{lbl}\n→ {wf}',
             ha='center', va='center', fontsize=7.5,
             color='#633806',
             bbox=dict(boxstyle='round,pad=0.2', fc='#FAEEDA',
                       ec='#854F0B', alpha=0.85))

# title and legend
ax1.set_title('Real space square lattice: nearest (NN) and next-nearest (NNN) neighbors\n'
              r'$k = 2\pi/a$,  all 8 directions and their wave function contributions',
              fontsize=11, pad=12)

legend_handles = [
    mpatches.Patch(fc='#534AB7', ec='#3C3489', label='Central atom (origin)'),
    mpatches.Patch(fc='#E6F1FB', ec='#185FA5',
                   label=r'NN (j=1):  distance $a$,  solid arrows'),
    mpatches.Patch(fc='#FAEEDA', ec='#854F0B',
                   label=r'NNN (j=2):  distance $a\sqrt{2}$,  dashed arrows'),
]
ax1.legend(handles=legend_handles, loc='upper left', fontsize=9,
           framealpha=0.9, bbox_to_anchor=(-0.02, 1.01))

plt.tight_layout()
fig1.savefig(os.path.join(SAVE_DIR, 'fig1_real_space_labeled.png'),
             dpi=200, bbox_inches='tight')
print("Saved fig1")

# ═══════════════════════════════════════════════════════════════════════════════
# Figure 2: j=1 mode density — annotated
# ═══════════════════════════════════════════════════════════════════════════════
fig2, ax2 = plt.subplots(figsize=(6, 6))
im2 = ax2.imshow(J1, extent=[0, 2*a, 0, 2*a], origin='lower',
                 cmap='RdBu_r', vmin=-4, vmax=4, aspect='equal')
add_colorbar(fig2, ax2, im2, 'density (a.u.)')

for (sx, sy) in ATOM_SITES:
    if sx <= 2*a and sy <= 2*a:
        ax2.plot(sx, sy, '^', ms=11, color='white',
                 mec='black', mew=0.8, zorder=5)
        ax2.text(sx+0.04*a, sy+0.07*a, 'atom\n+4',
                 fontsize=7, color='white', fontweight='bold',
                 bbox=dict(boxstyle='round,pad=0.15', fc='steelblue',
                           alpha=0.75, ec='none'))

ec = [(a/2,0),(0,a/2),(a,a/2),(a/2,a),(3*a/2,0),(2*a,a/2),
      (a/2,2*a),(3*a/2,2*a),(2*a,3*a/2),(0,3*a/2)]
for (sx, sy) in ec:
    if 0<=sx<=2*a and 0<=sy<=2*a:
        val = float(2*np.cos(k*sx) + 2*np.cos(k*sy))
        ax2.plot(sx, sy, 'x', ms=7, color='cyan', mew=1.8, zorder=5)
        ax2.text(sx+0.04*a, sy+0.04*a, f'{val:.0f}',
                 fontsize=7, color='cyan')

ax2.set_xlabel('x', fontsize=11); ax2.set_ylabel('y', fontsize=11)
ax2.set_title('j=1 mode:  $2\\cos(kx)+2\\cos(ky)$\n'
              r'▲ = atom site (+4)   × = edge center (0)', fontsize=10)
ax2.set_xticks([0,a,2*a]); ax2.set_xticklabels(['0','a','2a'])
ax2.set_yticks([0,a,2*a]); ax2.set_yticklabels(['0','a','2a'])
plt.tight_layout()
fig2.savefig(os.path.join(SAVE_DIR, 'fig2_mode_j1_annotated.png'), dpi=200)
print("Saved fig2")

# ═══════════════════════════════════════════════════════════════════════════════
# Figure 3: j=2 mode density — annotated
# ═══════════════════════════════════════════════════════════════════════════════
fig3, ax3 = plt.subplots(figsize=(6, 6))
im3 = ax3.imshow(J2, extent=[0, 2*a, 0, 2*a], origin='lower',
                 cmap='RdBu_r', vmin=-4, vmax=4, aspect='equal')
add_colorbar(fig3, ax3, im3, 'density (a.u.)')

for (sx, sy) in ATOM_SITES:
    if sx <= 2*a and sy <= 2*a:
        ax3.plot(sx, sy, '^', ms=11, color='white',
                 mec='black', mew=0.8, zorder=5)
        ax3.text(sx+0.04*a, sy+0.07*a, 'atom\n+4',
                 fontsize=7, color='white', fontweight='bold',
                 bbox=dict(boxstyle='round,pad=0.15', fc='steelblue',
                           alpha=0.75, ec='none'))

for (sx, sy) in FACE_SITES:
    ax3.plot(sx, sy, '*', ms=10, color='orange',
             mec='white', mew=0.5, zorder=5)
    ax3.text(sx+0.04*a, sy-0.16*a, 'face\n+4',
             fontsize=7, color='orange',
             bbox=dict(boxstyle='round,pad=0.1', fc='black',
                       alpha=0.45, ec='none'))

for (sx, sy) in ec:
    if 0<=sx<=2*a and 0<=sy<=2*a:
        val = float(4*np.cos(k*sx)*np.cos(k*sy))
        ax3.plot(sx, sy, 'x', ms=7, color='cyan', mew=1.8, zorder=5)
        ax3.text(sx+0.04*a, sy+0.04*a, f'{val:.0f}',
                 fontsize=7, color='cyan')

ax3.set_xlabel('x', fontsize=11); ax3.set_ylabel('y', fontsize=11)
ax3.set_title('j=2 mode:  $4\\cos(kx)\\cdot\\cos(ky)$\n'
              r'▲ = atom site (+4)   ★ = face center (+4)   × = edge center (−4)',
              fontsize=10)
ax3.set_xticks([0,a,2*a]); ax3.set_xticklabels(['0','a','2a'])
ax3.set_yticks([0,a,2*a]); ax3.set_yticklabels(['0','a','2a'])
plt.tight_layout()
fig3.savefig(os.path.join(SAVE_DIR, 'fig3_mode_j2_annotated.png'), dpi=200)
print("Saved fig3")

# ═══════════════════════════════════════════════════════════════════════════════
# Figure 4: Combined density — three values of A2
# ═══════════════════════════════════════════════════════════════════════════════
fig4, axes = plt.subplots(1, 3, figsize=(15, 5))
for idx, (A1v, A2v) in enumerate([(0.5,0.0),(0.5,0.2),(0.5,0.4)]):
    ax = axes[idx]
    Z  = density(X, Y, k, A1v, A2v)
    im = ax.imshow(Z, extent=[0,2*a,0,2*a], origin='lower',
                   cmap='RdBu_r', vmin=Z.min(), vmax=Z.max(), aspect='equal')
    add_colorbar(fig4, ax, im, 'n(r)')

    for (sx, sy) in ATOM_SITES:
        if sx<=2*a and sy<=2*a:
            ax.plot(sx, sy, '^', ms=11, color='white',
                    mec='black', mew=0.8, zorder=6)
            val = density(np.array([[sx]]), np.array([[sy]]),
                          k, A1v, A2v)[0,0]
            ax.text(sx+0.05*a, sy+0.08*a,
                    f'atom\n{val:.2f}', fontsize=6.5, color='white',
                    fontweight='bold',
                    bbox=dict(boxstyle='round,pad=0.15', fc='#1a4a7a',
                              alpha=0.75, ec='none'))

    for (sx, sy) in FACE_SITES:
        ax.plot(sx, sy, 'D', ms=7, color='orange',
                mec='white', mew=0.6, zorder=5)
        val = density(np.array([[sx]]), np.array([[sy]]),
                      k, A1v, A2v)[0,0]
        ax.text(sx+0.04*a, sy-0.15*a,
                f'NNN\n{val:.2f}', fontsize=6.5, color='orange',
                bbox=dict(boxstyle='round,pad=0.15', fc='black',
                          alpha=0.45, ec='none'))

    ax.set_title(f'$A_1={A1v}$,  $A_2={A2v}$', fontsize=11)
    ax.set_xlabel('x', fontsize=10); ax.set_ylabel('y', fontsize=10)
    ax.set_xticks([0,a,2*a]); ax.set_xticklabels(['0','a','2a'])
    ax.set_yticks([0,a,2*a]); ax.set_yticklabels(['0','a','2a'])

legend_elements = [
    plt.Line2D([0],[0], marker='^', color='w', markerfacecolor='white',
               mec='black', ms=10, label='Atom (lattice site)'),
    plt.Line2D([0],[0], marker='D', color='w', markerfacecolor='orange',
               ms=8, label='NNN face center'),
]
fig4.legend(handles=legend_elements, loc='lower center', ncol=2,
            fontsize=9, bbox_to_anchor=(0.5,-0.04))
fig4.suptitle(
    r'$n(\mathbf{r}) = A_1[2\cos(kx)+2\cos(ky)] + A_2\cdot4\cos(kx)\cos(ky)$'
    '\n  increasing $A_2$ sharpens atomic peaks',
    fontsize=11, y=1.02)
plt.tight_layout()
fig4.savefig(os.path.join(SAVE_DIR,'fig4_combined_annotated.png'),
             dpi=200, bbox_inches='tight')
print("Saved fig4")

# ═══════════════════════════════════════════════════════════════════════════════
# Figure 5: 1D cross-section
# ═══════════════════════════════════════════════════════════════════════════════
fig5, ax5 = plt.subplots(figsize=(9, 4))
x1d = np.linspace(0, 2*a, 1000)
j1_1d  = mode_j1(x1d, np.zeros_like(x1d), k)
j2_1d  = mode_j2(x1d, np.zeros_like(x1d), k)
sum_1d = A1*j1_1d + A2*j2_1d

ax5.plot(x1d/a, j1_1d,  '-',  color='steelblue',  lw=2,   label=f'j=1 (A₁={A1})')
ax5.plot(x1d/a, j2_1d,  '--', color='darkorange', lw=2,   label=f'j=2 (A₂={A2})')
ax5.plot(x1d/a, sum_1d, '-',  color='purple',     lw=2.5, label='combined n(x,0)')

for ix in range(3):
    ax5.axvline(ix, color='gray', lw=0.8, ls=':', alpha=0.6)
    yv = density(np.array([[ix*a]]), np.zeros((1,1)), k, A1, A2)[0,0]
    ax5.plot(ix, yv, '^', ms=10, color='purple',
             mec='black', mew=0.8, zorder=6)
    ax5.text(ix+0.03, yv+0.05, 'atom', fontsize=8, color='purple')

for ix in [0.5, 1.5]:
    yv = density(np.array([[ix*a]]), np.zeros((1,1)), k, A1, A2)[0,0]
    ax5.plot(ix, yv, 'x', ms=9, color='darkorange', mew=2, zorder=6)
    ax5.text(ix+0.03, yv-0.18, f'{yv:.2f}', fontsize=8, color='darkorange')

ax5.axhline(0, color='gray', lw=0.5)
ax5.set_xlabel('x / a', fontsize=11)
ax5.set_ylabel('n(x, y=0)', fontsize=11)
ax5.set_title('1D cross-section y=0: j=1 vs j=2 contributions', fontsize=10)
ax5.legend(fontsize=9)
ax5.grid(True, alpha=0.2)
plt.tight_layout()
fig5.savefig(os.path.join(SAVE_DIR,'fig5_1d_crosssection.png'), dpi=200)
print("Saved fig5")

plt.show()
print("\nAll figures saved to:", SAVE_DIR)
