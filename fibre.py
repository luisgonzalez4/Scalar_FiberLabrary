import numpy as np
import matplotlib.pyplot as plt
from scipy.special import j0

# Parámetro normalizado típico del modo fundamental (primer cero de J0)
pa = 2.405  

# Eje normalizado rho/a
x = np.linspace(0, 3, 500)

# Funciones
J = j0(pa * x)                 # J0(p*rho) con rho = a*x
w_over_a = 0.65
G = np.exp(-(x / w_over_a)**2)  # Gaussiana aproximada

# Gráfica
plt.figure(figsize=(5,4))
plt.plot(x, J, label=r"$J_0(p\rho)$")
plt.plot(x, G, label="Gaussiana")
plt.xlabel(r"$\rho / a$")
plt.ylabel("Amplitud (normalizada)")
plt.legend()
plt.tight_layout()
plt.show()
