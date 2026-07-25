#!/usr/bin/env python3
"""Visualize Z* cubic residual f*(Z*) = Z*³ + d2·Z*² + d1·Z* + d0
across a range of pressures and compositions to build intuition
about the root structure in Z* = Z - B space.

For PR EOS: d2 = 4B-1, d1 = A+2B(B-2), d0 = -2B²
Guaranteed: f*(0) = -2B² < 0, f*(1) = A >= 0
So all physical roots lie in (0, 1].
"""

import numpy as np
import matplotlib.pyplot as plt
import sys, os
sys.path.insert(0, os.path.dirname(__file__))

from pyrestoolbox.gas.gas import (
    gas_tc_pc, _BNS_TCS, _BNS_PCS, _BNS_ACF,
    _BNS_OMEGAA, _BNS_OMEGAB, _calc_bips_fast,
)

DEGF2R = 459.67


def get_AB(p_psia, degf, sg, co2=0, h2s=0, n2=0, h2=0):
    """Compute EOS mixing parameters A and B for given conditions."""
    tc, pc = gas_tc_pc(sg, co2, h2s, n2, h2, 'BNS')
    degR = degf + DEGF2R

    tcs = _BNS_TCS.copy()
    pcs = _BNS_PCS.copy()
    tcs[-1], pcs[-1] = tc, pc

    z = np.array([co2, h2s, n2, h2, 1 - co2 - h2s - n2 - h2])
    trs = degR / tcs
    m = 0.37464 + 1.54226 * _BNS_ACF - 0.26992 * _BNS_ACF**2
    alpha = (1 + m * (1 - np.sqrt(trs)))**2
    kij = _calc_bips_fast(degR, tc)

    prs = p_psia / pcs
    Ai = _BNS_OMEGAA * alpha * prs / trs**2
    Bi = _BNS_OMEGAB * prs / trs

    sqrt_Ai = np.sqrt(Ai)
    w = z * sqrt_Ai
    onemk = 1.0 - kij
    A = float(np.sum((w @ onemk) * w))
    B = float(Bi @ z)
    return A, B


def zstar_cubic(zs, A, B):
    """Evaluate f*(Z*) = Z*³ + d2·Z*² + d1·Z* + d0"""
    d2 = 4 * B - 1
    d1 = A + 2 * B * (B - 2)
    d0 = -2 * B * B
    return zs**3 + d2 * zs**2 + d1 * zs + d0


def zstar_deriv(zs, A, B):
    """Evaluate f*'(Z*) = 3Z*² + 2d2·Z* + d1"""
    d2 = 4 * B - 1
    d1 = A + 2 * B * (B - 2)
    return 3 * zs**2 + 2 * d2 * zs + d1


zs_range = np.linspace(0, 1, 500)

# --- Plot 1: Varying pressure (sweet gas) ---
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

pressures = [500, 2000, 5000, 10000, 15000, 20000, 30000]
colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(pressures)))

ax = axes[0, 0]
for p, c in zip(pressures, colors):
    A, B = get_AB(p, 200, 0.7)
    f = zstar_cubic(zs_range, A, B)
    ax.plot(zs_range, f, color=c, label=f'{p} psia (A={A:.2f}, B={B:.2f})')
ax.axhline(0, color='k', linewidth=0.5, linestyle='--')
ax.set_xlabel('Z*')
ax.set_ylabel('f*(Z*)')
ax.set_title('SG=0.7, T=200°F — Varying Pressure')
ax.legend(fontsize=7, loc='upper left')
ax.set_xlim(0, 1)
ax.grid(True, alpha=0.3)

# --- Plot 2: Varying temperature ---
ax = axes[0, 1]
temperatures = [80, 150, 200, 300, 400]
colors_t = plt.cm.coolwarm(np.linspace(0.1, 0.9, len(temperatures)))
for t, c in zip(temperatures, colors_t):
    A, B = get_AB(5000, t, 0.7)
    f = zstar_cubic(zs_range, A, B)
    ax.plot(zs_range, f, color=c, label=f'T={t}°F (A={A:.2f}, B={B:.2f})')
ax.axhline(0, color='k', linewidth=0.5, linestyle='--')
ax.set_xlabel('Z*')
ax.set_ylabel('f*(Z*)')
ax.set_title('SG=0.7, P=5000 psia — Varying Temperature')
ax.legend(fontsize=7, loc='upper left')
ax.set_xlim(0, 1)
ax.grid(True, alpha=0.3)

# --- Plot 3: Varying composition (impurities at 10000 psia) ---
ax = axes[1, 0]
cases = [
    ('Pure HC', dict(co2=0, h2s=0, n2=0, h2=0)),
    ('20% CO₂', dict(co2=0.2, h2s=0, n2=0, h2=0)),
    ('10% H₂S', dict(co2=0, h2s=0.1, n2=0, h2=0)),
    ('10% N₂', dict(co2=0, h2s=0, n2=0.1, h2=0)),
    ('30% H₂', dict(co2=0, h2s=0, n2=0, h2=0.3)),
    ('20% CO₂ + 10% H₂S', dict(co2=0.2, h2s=0.1, n2=0, h2=0)),
]
colors_c = plt.cm.Set1(np.linspace(0, 0.8, len(cases)))
for (label, kwargs), c in zip(cases, colors_c):
    A, B = get_AB(10000, 200, 0.7, **kwargs)
    f = zstar_cubic(zs_range, A, B)
    ax.plot(zs_range, f, color=c, label=f'{label} (A={A:.2f}, B={B:.2f})')
ax.axhline(0, color='k', linewidth=0.5, linestyle='--')
ax.set_xlabel('Z*')
ax.set_ylabel('f*(Z*)')
ax.set_title('SG=0.7, P=10000 psia, T=200°F — Varying Composition')
ax.legend(fontsize=7, loc='upper left')
ax.set_xlim(0, 1)
ax.grid(True, alpha=0.3)

# --- Plot 4: f* and f*' together for a near-critical case ---
ax = axes[1, 1]
# 80% CO2 at 90°F — CO2 Tc≈88°F so Tr≈1.0, genuinely near-critical
# 80% CO2 + 10% H2S at 40°F — genuinely near-critical for CO2 (Tc≈88°F)
# with well-separated 3-root VLE behavior in Z* space
p_nc = 600
A, B = get_AB(p_nc, 40, 0.7, co2=0.8, h2s=0.1)

f = zstar_cubic(zs_range, A, B)
fp = zstar_deriv(zs_range, A, B)

ax.plot(zs_range, f, 'b-', linewidth=2, label="f*(Z*)")
ax.plot(zs_range, fp, 'r--', linewidth=1.5, label="f*'(Z*)")
ax.axhline(0, color='k', linewidth=0.5, linestyle='--')

# Mark roots
roots = []
for i in range(len(zs_range) - 1):
    if f[i] * f[i + 1] < 0:
        root = zs_range[i] - f[i] * (zs_range[i + 1] - zs_range[i]) / (f[i + 1] - f[i])
        roots.append(root)
for r in roots:
    ax.axvline(r, color='green', linewidth=0.8, linestyle=':', alpha=0.7)
    ax.plot(r, 0, 'go', markersize=8)

ax.set_xlabel('Z*')
ax.set_ylabel('f*(Z*) / f*\'(Z*)')
ax.set_title(f'Near-critical: 80% CO₂ + 10% H₂S, P={p_nc} psia, T=40°F\n'
             f'A={A:.2f}, B={B:.2f} — '
             f'{len(roots)} root(s) in (0,1]: {[f"{r:.4f}" for r in roots]}')
ax.legend(fontsize=9)
ax.set_xlim(0, 1)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('/home/mark/projects/pyResToolbox/zstar_residual_overview.png', dpi=150)
plt.close()
print('Saved zstar_residual_overview.png')

# --- Plot 5: Newton convergence from Z*=1 ---
fig, axes = plt.subplots(1, 3, figsize=(16, 5))

test_cases = [
    ('Low P: 500 psia, 200°F', 500, 200, dict()),
    ('High P: 20000 psia, 200°F', 20000, 200, dict()),
    ('Near-critical: 80% CO₂ + 10% H₂S\n600 psia, 40°F', 600, 40, dict(co2=0.8, h2s=0.1)),
]

for ax, (title, p, t, comp) in zip(axes, test_cases):
    A, B = get_AB(p, t, 0.7, **comp)

    # Show the cubic
    f = zstar_cubic(zs_range, A, B)
    ax.plot(zs_range, f, 'b-', linewidth=2, label='f*(Z*)')
    ax.axhline(0, color='k', linewidth=0.5, linestyle='--')

    # Trace Newton iterations from Z*=1
    zs = 1.0
    d2 = 4 * B - 1
    d1 = A + 2 * B * (B - 2)
    d0 = -2 * B * B

    iterates = [zs]
    for _ in range(15):
        fv = zs**3 + d2 * zs**2 + d1 * zs + d0
        fpv = 3 * zs**2 + 2 * d2 * zs + d1
        if abs(fpv) < 1e-30:
            break
        dz = fv / fpv
        zs -= dz
        iterates.append(zs)
        if abs(dz) < 1e-12:
            break

    # Plot iteration path
    for i, zsi in enumerate(iterates):
        fv = zsi**3 + d2 * zsi**2 + d1 * zsi + d0
        ax.plot(zsi, fv, 'ro' if i == 0 else 'g^' if i == len(iterates) - 1 else 'k.',
                markersize=10 if i in (0, len(iterates) - 1) else 5, zorder=5)
        if i < len(iterates) - 1:
            # Draw Newton tangent line to x-axis
            fpv = 3 * zsi**2 + 2 * d2 * zsi + d1
            next_zs = zsi - fv / fpv
            ax.plot([zsi, next_zs], [fv, 0], 'k-', linewidth=0.5, alpha=0.5)
            ax.plot([next_zs, next_zs], [0, next_zs**3 + d2 * next_zs**2 + d1 * next_zs + d0],
                    'k:', linewidth=0.5, alpha=0.5)

    ax.set_xlabel('Z*')
    ax.set_ylabel('f*(Z*)')
    ax.set_title(f'{title}\nA={A:.3f}, B={B:.4f}, {len(iterates)-1} iterations')
    ax.set_xlim(0, 1.05)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # Annotate start and end
    ax.annotate(f'Start: Z*=1', xy=(iterates[0], zstar_cubic(iterates[0], A, B)),
                xytext=(0.5, 0.85), textcoords='axes fraction',
                arrowprops=dict(arrowstyle='->', color='red'), fontsize=8, color='red')
    ax.annotate(f'Root: Z*={iterates[-1]:.6f}', xy=(iterates[-1], 0),
                xytext=(0.5, 0.15), textcoords='axes fraction',
                arrowprops=dict(arrowstyle='->', color='green'), fontsize=8, color='green')

plt.tight_layout()
plt.savefig('/home/mark/projects/pyResToolbox/zstar_newton_convergence.png', dpi=150)
plt.close()
print('Saved zstar_newton_convergence.png')

# --- Plot 6: Comparison of Z-space vs Z*-space cubics at high pressure ---
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

p, t = 20000, 200
A, B = get_AB(p, t, 0.7)

# Z-space cubic
c2 = -(1 - B)
c1 = A - 3 * B**2 - 2 * B
c0 = -(A * B - B**2 - B**3)

z_range = np.linspace(-3, 4, 500)
fz = z_range**3 + c2 * z_range**2 + c1 * z_range + c0

ax = axes[0]
ax.plot(z_range, fz, 'b-', linewidth=2)
ax.axhline(0, color='k', linewidth=0.5, linestyle='--')
ax.axvline(B, color='orange', linewidth=1, linestyle='--', label=f'B = {B:.3f}')

# Mark inflection + 1 (old starting point)
z0_old = -c2 / 3 + 1
ax.plot(z0_old, z0_old**3 + c2 * z0_old**2 + c1 * z0_old + c0, 'rs',
        markersize=12, label=f'Old start: -c2/3+1 = {z0_old:.2f}', zorder=5)

# Mark roots
for i in range(len(z_range) - 1):
    if fz[i] * fz[i + 1] < 0:
        root = z_range[i] - fz[i] * (z_range[i + 1] - z_range[i]) / (fz[i + 1] - fz[i])
        ax.plot(root, 0, 'go', markersize=8, zorder=5)
        ax.annotate(f'Z={root:.2f}', xy=(root, 0), xytext=(root, 5),
                    fontsize=8, ha='center')

ax.set_xlabel('Z')
ax.set_ylabel('f(Z)')
ax.set_title(f'Z-space cubic at {p} psia\nc2={c2:.3f}, c1={c1:.3f}, c0={c0:.3f}')
ax.legend(fontsize=8)
ax.set_ylim(-30, 30)
ax.grid(True, alpha=0.3)

# Z*-space cubic
ax = axes[1]
f_zs = zstar_cubic(zs_range, A, B)
ax.plot(zs_range, f_zs, 'b-', linewidth=2)
ax.axhline(0, color='k', linewidth=0.5, linestyle='--')
ax.plot(1.0, zstar_cubic(1.0, A, B), 'gs', markersize=12,
        label=f'Start: Z*=1, f*={A:.3f}', zorder=5)

for i in range(len(zs_range) - 1):
    if f_zs[i] * f_zs[i + 1] < 0:
        root = zs_range[i] - f_zs[i] * (zs_range[i + 1] - zs_range[i]) / (f_zs[i + 1] - f_zs[i])
        ax.plot(root, 0, 'go', markersize=8, zorder=5)
        ax.annotate(f'Z*={root:.4f}', xy=(root, 0), xytext=(root + 0.05, -0.5),
                    fontsize=8)

ax.set_xlabel('Z*')
ax.set_ylabel('f*(Z*)')
d2 = 4 * B - 1
d1_val = A + 2 * B * (B - 2)
d0_val = -2 * B * B
ax.set_title(f'Z*-space cubic at {p} psia\nd2={d2:.3f}, d1={d1_val:.3f}, d0={d0_val:.3f}')
ax.legend(fontsize=8)
ax.set_xlim(0, 1.05)
ax.grid(True, alpha=0.3)

plt.suptitle(f'Why Z* works: P={p} psia, T={t}°F, SG=0.7 (A={A:.3f}, B={B:.4f})',
             fontsize=13, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig('/home/mark/projects/pyResToolbox/zstar_vs_z_comparison.png', dpi=150, bbox_inches='tight')
plt.close()
print('Saved zstar_vs_z_comparison.png')
