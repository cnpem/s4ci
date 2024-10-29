
import numpy as np
import form_factors as f
from random import randint

rand = randint(1, 4)
rand=4

if rand == 1:
    ff = f.Esfera(50, 1)
if rand == 2:
    ff = f.Elipsoide(50, 1.5)
if rand == 3:
    ff = f.Cilindro(50, 200)
if rand == 4:
    ff = f.CascaEsferica(50, 75, 1.3, 1)
q = np.logspace(-2.4, 0, 500)
ff.setq(q)
ff.ff()
#ff.distribucao(5, 0, 3, 'ln')
#ff.intensidade = ff.ff_dist()
f.ruido(ff, 0, 1e6)

ff2 = f.Esfera(30, 0.75)
ff2.setq(q)
ff2.ff()
f.ruido(ff2, 0, 1e6)

#f = open("dado_simulado.txt", "w")
#f.close()
data = np.column_stack([ff.q, ff.intensidade+ff2.intensidade])
np.savetxt(r'C:\Users\pedro.guidolim\OneDrive - CNPEM - Centro Nacional de Pesquisa em Energia e Materiais\Documentos\dado_simulado.txt', data)
print(rand)