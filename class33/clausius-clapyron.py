import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

# 1. Generate Data: Temperature range from 25 C to 100 C
T_C = np.linspace(25, 100, 20)
T_K = T_C + 273.15

# Antoine equation constants for water (T in Celsius, P in mmHg)
A = 8.07131
B = 1730.63
C = 233.426

# Calculate Saturation Pressure in mmHg, then convert to Pa
P_mmHg = 10**(A - B / (T_C + C))
P_Pa = P_mmHg * 133.322

# 2. Prepare Clausius-Clapeyron variables
x = 1 / T_K
y = np.log(P_Pa)

# 3. Perform Linear Regression
slope, intercept, r_value, p_value, std_err = linregress(x, y)

# 4. Extract Enthalpy of Vaporization
# The Clausius-Clapeyron equation dictates: slope = - deltaH_vap / R
R = 8.314 # Universal Gas Constant in J/(mol K)
deltaH_vap_J_mol = -slope * R
deltaH_vap_kJ_mol = deltaH_vap_J_mol / 1000

print(f"Calculated Enthalpy of Vaporization: {deltaH_vap_kJ_mol:.2f} kJ/mol")

# 5. Create the Plot
plt.figure(figsize=(8, 6))
plt.scatter(x, y, color='blue', label='Data', zorder=5)

# Plot the line of best fit
plt.plot(x, slope * x + intercept, color='red', linestyle='--', 
         label=f'Linear Fit: slope = {slope:.1f} K')

# Formatting for textbook/lecture quality
plt.title(r'Clausius-Clapeyron Plot for Water', fontsize=16)
plt.xlabel(r'$1/T$ (K$^{-1}$)', fontsize=14)
plt.ylabel(r'$\ln P^{sat}$ ($P$ in Pa)', fontsize=14)
plt.grid(True, linestyle=':', alpha=0.7)
plt.legend(fontsize=12)

# Annotate the calculated delta H on the plot
annotation_text = (r'Slope $= -\frac{\Delta H_{vap}}{R}$' + '\n' +
                   r'$\Delta H_{vap} \approx$ ' + f'{deltaH_vap_kJ_mol:.1f} kJ/mol')
plt.text(0.05, 0.15, annotation_text, transform=plt.gca().transAxes,
         fontsize=12, bbox=dict(facecolor='white', alpha=0.8, edgecolor='gray'))

plt.tight_layout()
plt.savefig('clausius_clapeyron_water.png', dpi=300)
plt.show()