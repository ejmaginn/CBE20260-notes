import numpy as np
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(5.5, 5))

# ── Vapor-liquid dome ──────────────────────────────────────────────────────────
Vc, Pc = 0.50, 1.00
t = np.linspace(0, 1, 200)

V_vap_curve = Vc + 1.60 * t**0.65    # right (vapor) saturation branch
V_liq_curve = Vc - 0.32 * t**0.55    # left  (liquid) saturation branch
P_dome      = Pc * (1 - t)

ax.plot(V_liq_curve, P_dome, 'k-', lw=1.8)
ax.plot(V_vap_curve, P_dome, 'k-', lw=1.8)
ax.plot(Vc, Pc, 'ko', ms=5)

# ── Key pressures & saturation volumes ────────────────────────────────────────
P_sat = 0.38
P_D   = 0.80

t_sat = 1 - P_sat / Pc
V_B   = Vc + 1.60 * t_sat**0.65   # B = saturated VAPOR  (right branch)
V_C   = Vc - 0.32 * t_sat**0.55   # C = saturated LIQUID (left branch)
V_D   = V_C - 0.03                 # D = compressed liquid (nearly same V as C)

# A = low-pressure vapor, well to the right
P_A = 0.06
V_A = V_B + 1.20

# ── Leg 1: vapor isotherm from A → B ─────────────────────────────────────────
# Use power-law P = P_A * (V_A / V)^n, calibrated to pass through B at P_sat
n_exp = np.log(P_sat / P_A) / np.log(V_A / V_B)
V_leg1 = np.linspace(V_A, V_B, 200)
P_leg1 = P_A * (V_A / V_leg1) ** n_exp

ax.plot(V_leg1, P_leg1, color='royalblue', lw=2.0, zorder=4)

# Arrow on Leg 1 (near B)
idx = 170
ax.annotate('', xy=(V_leg1[idx+5], P_leg1[idx+5]),
            xytext=(V_leg1[idx], P_leg1[idx]),
            arrowprops=dict(arrowstyle='->', color='royalblue', lw=1.8))

# ── Leg 2: tie line B → C (right to left at P_sat) ───────────────────────────
ax.plot([V_B, V_C], [P_sat, P_sat], color='royalblue', lw=2.0, zorder=4)
ax.annotate('', xy=(V_C + 0.05, P_sat),
            xytext=(V_C + 0.25, P_sat),
            arrowprops=dict(arrowstyle='->', color='royalblue', lw=1.8))

# ── Leg 3: liquid from C → D (nearly vertical, steep) ────────────────────────
ax.plot([V_C, V_D], [P_sat, P_D], color='royalblue', lw=2.0, zorder=4)
ax.annotate('', xy=(V_D, P_D - 0.04),
            xytext=(V_D, P_D - 0.18),
            arrowprops=dict(arrowstyle='->', color='royalblue', lw=1.8))

# ── Point markers ─────────────────────────────────────────────────────────────
mkw = dict(ms=7, zorder=6, clip_on=False, color='royalblue')
for pt in [(V_A, P_A), (V_B, P_sat), (V_C, P_sat), (V_D, P_D)]:
    ax.plot(*pt, 'o', **mkw)

# ── Point labels ──────────────────────────────────────────────────────────────
ax.text(V_A + 0.05, P_A,        r'$A\ (P \approx 0)$',
        va='center', fontsize=10)
ax.text(V_B + 0.04, P_sat - 0.05, r'$B$',
        va='top', ha='left', fontsize=11)
ax.text(V_C - 0.04, P_sat - 0.05, r'$C$',
        va='top', ha='right', fontsize=11)
ax.text(V_D - 0.04, P_D,         r'$D$',
        va='center', ha='right', fontsize=11)

# ── Leg labels ────────────────────────────────────────────────────────────────
mid1_idx = 110
ax.text(V_leg1[mid1_idx] + 0.08, P_leg1[mid1_idx],
        'Leg 1\n(vapor)', color='royalblue', fontsize=8.5,
        ha='left', va='center')

ax.text((V_B + V_C) / 2, P_sat + 0.05,
        r'Leg 2  (VLE, $f^V\!=\!f^L$)', color='royalblue',
        fontsize=8.5, ha='center', va='bottom')

ax.text(V_D - 0.08, (P_sat + P_D) / 2,
        'Leg 3\n(liquid)', color='royalblue', fontsize=8.5,
        ha='right', va='center')

# ── Pressure dashed guidelines ────────────────────────────────────────────────
ax.axhline(P_sat, ls='--', color='gray', lw=0.8, zorder=1)
ax.axhline(P_D,   ls='--', color='gray', lw=0.8, zorder=1)
ax.set_yticks([P_sat, P_D])
ax.set_yticklabels([r'$P^{\mathrm{sat}}$', r'$P$'], fontsize=11)

# ── Phase region labels ───────────────────────────────────────────────────────
ax.text(0.13, 0.55, 'Liquid', fontsize=9, color='0.45',
        transform=ax.transAxes, ha='center')
ax.text(0.80, 0.55, 'Vapor', fontsize=9, color='0.45',
        transform=ax.transAxes, ha='center')
ax.text(0.47, 0.30, 'Two-phase', fontsize=9, color='0.45',
        transform=ax.transAxes, ha='center', style='italic')
ax.text(Vc + 0.05, Pc + 0.04, 'Critical\npoint', fontsize=8, va='bottom')

# ── Axes cosmetics ────────────────────────────────────────────────────────────
ax.set_xlim(0.05, V_A + 0.30)
ax.set_ylim(-0.02, 1.12)
ax.set_xlabel(r'Molar volume, $V$', fontsize=12)
ax.set_ylabel(r'Pressure, $P$', fontsize=12)
ax.set_xticks([])
ax.spines[['top', 'right']].set_visible(False)

fig.tight_layout()
fig.savefig(
    '/Users/edwardmaginn/Thermo-2026/Lecture-Notes-2026/class34/images/liquid-fugacity-path.png',
    dpi=180, bbox_inches='tight'
)
print('saved.')
