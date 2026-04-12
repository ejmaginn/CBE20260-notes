import numpy as np
import matplotlib.pyplot as plt

# Simulate some realistic data for (Z-1)*V vs 1/V
rho = np.linspace(0, 10, 100) # Density (1/V) in mol/L

# T1: Low Temperature (Negative B)
B1 = -0.15
C1 = 0.005
y1 = B1 + C1 * rho

# T2: Boyle Temperature (B = 0)
B2 = 0.0
C2 = 0.003
y2 = B2 + C2 * rho

# T3: High Temperature (Positive B)
B3 = 0.05
C3 = 0.001
y3 = B3 + C3 * rho

plt.figure(figsize=(8, 6))
plt.plot(rho, y1, 'b-', linewidth=2, label=r'Low $T$ ($T < T_{Boyle}$)')
plt.plot(rho, y2, 'g-', linewidth=2, label=r'Boyle Temp ($T = T_{Boyle}$)')
plt.plot(rho, y3, 'r-', linewidth=2, label=r'High $T$ ($T > T_{Boyle}$)')

# Formatting to look like a textbook
plt.axhline(0, color='black', linewidth=1)
plt.axvline(0, color='black', linewidth=1)
plt.xlim(0, 10)
plt.ylim(-0.2, 0.30)

# Annotate Intercepts (B)
plt.plot(0, B1, 'bo', markersize=8)
plt.plot(0, B2, 'go', markersize=8)
plt.plot(0, B3, 'ro', markersize=8)

plt.annotate(r'$B(T_{low})$', xy=(0.2, B1-0.01), color='blue', fontsize=12)
plt.annotate(r'$B(T_{high})$', xy=(0.2, B3+0.01), color='red', fontsize=12)

# Removed the \underline command here!
plt.xlabel(r'Molar Density, $\rho = 1/V$', fontsize=20)
plt.ylabel(r'$(Z-1)V$', fontsize=20)
plt.title('Experimental Determination of Virial Coefficients', fontsize=20)
plt.legend(fontsize=20)
plt.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()
plt.show()