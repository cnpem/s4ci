import numpy as np
import form_factors as f
import matplotlib.pyplot as plt

q = np.logspace(-3, 0, 100)

form_factor = f.Esfera(50, 1.2)
form_factor.setq(q)
intensidade = form_factor.ff()

fig, ax = plt.subplots()
ax.plot(q, intensidade)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylabel('Intensity [u.a.]')
ax.set_xlabel('Q '+'$[\u212B^{-1}]$') #u212B == Angstrom
plt.show()