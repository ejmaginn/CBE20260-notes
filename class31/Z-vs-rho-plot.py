import numpy as np
import matplotlib.pyplot as plt

# Properties for a representative fluid: Methane
Tc = 190.6  # K
Pc = 45.99e5 # Pa
omega = 0.011
R = 8.314 # J/(mol K)

# Peng-Robinson parameters
a = 0.45724 * (R * Tc)**2 / Pc
b = 0.07780 * R * Tc / Pc
kappa = 0.37464 + 1.54226 * omega - 0.26992 * omega**2

def calc_Z(T, rho_mol_L):
    # Convert density from mol/L to mol/m^3, then volume to m^3/mol
    v = 1.0 / (rho_mol_L * 1000.0) 
    alpha = (1 + kappa * (1 - np.sqrt(T / Tc)))**2
    
    # Peng-Robinson Equation
    P = (R * T) / (v - b) - (a * alpha) / (v**2 + 2*b*v - b**2)
    
    # Compressibility factor
    Z = P * v / (R * T)
    return Z

# Density range (gas to moderate density)
# 0 to 12 mol/L is a good range to show the curvature
rho = np.linspace(0.001, 12, 500) 

plt.figure(figsize=(8, 6))

# Plot supercritical/critical isotherms
T_r_values = [1.0, 1.2, 1.5, 2.0, 3.0]
colors = ['purple', 'blue', 'green', 'orange', 'red']

for Tr, color in zip(T_r_values, colors):
    T = Tr * Tc
    Z = calc_Z(T, rho)
    plt.plot(rho, Z, color=color, lw=2.5, label=f'$T_r = {Tr}$ ($T = {T:.1f}$ K)')

# Formatting to look like a textbook
plt.axhline(1.0, color='black', linewidth=1.5, linestyle='--')
plt.axvline(0, color='black', linewidth=1)
plt.xlim(0, 12)
plt.ylim(0.5, 1.4)

plt.xlabel(r'Molar Density, $\rho = 1/V$ (mol/L)', fontsize=14)
plt.ylabel('Compressibility Factor, $Z$', fontsize=14)
plt.title('Peng-Robinson $Z$ vs. Density (Methane)', fontsize=16)
plt.legend(fontsize=12, loc='lower right')
plt.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()

# Save the figure
plt.savefig('Z_vs_rho_PR.png', dpi=300)
plt.show()