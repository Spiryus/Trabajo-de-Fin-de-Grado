# -*- coding: utf-8 -*-


"""
Modelo numérico completo de interior estelar

Ajuste del radio y luminosidad totales del modelo numérico de interior estelar
para conseguir el menor error relativo total en la frontera radiativo-
convectiva.

Posteriormente, se muestran los parámetros de la estrella, los valores de los
mismos en cada capa, una gráfica con el ajuste del error y otra con los
distintos parámetros representados en función del radio.
"""


import numpy as np
from estrella import star
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from matplotlib.patches import ConnectionPatch


# ---------------------- Valores del modelo de prueba ---------------------- #
#M_tot = 5.0
#X = 0.75
#Y = 0.22
#R_tot = 11.5
#L_tot = 70.0
#T_c = 2.0

# --------------- Valores del modelo a resolver para el TFG ---------------- #
#M_tot = 4.8
#X = 0.75
#Y = 0.20
#R_tot = 12.0
#L_tot = 35.0
#T_c = 1.5
# ----------------------------------------------------------------------------

num_capas = 100

M_tot = float(input('Masa total (x 10^33 g): '))
X = float(input('Fracción en masa de hidrógeno: '))
Y = float(input('Fracción en masa de helio: '))
R_tot = float(input('Radio total (x 10^10 cm): '))
L_tot = float(input('Luminosidad total (x 10^33 erg/s): '))
T_c = float(input('Temperatura central (x 10^7 K): '))
num_capas = float(input('Número de capas (100 por defecto): '))

print('Calculando...')
print('\n')

# =============================================================================
#      AJUSTE DE RADIO Y LUMINOSIDAD TOTALES, R_tot y L_tot
# =============================================================================
dR = 0.5
dL = 5
R_init = R_tot - 3*dR
R_fin = R_tot + 3*dR
L_init = L_tot - 3*dL
L_fin = L_tot + 3*dL

i = 0
E_rT = np.zeros((7, 7))
for R in np.arange(R_init, R_fin+.1, dR):
    j = 0
    for L in np.arange(L_init, L_fin+.1, dL):
        param, E_rT_i, T_c_def = star(M_tot, X, Y, R, L, T_c, num_capas)
        E_rT[i, j] = E_rT_i
        j += 1
    i += 1
minimo = np.min(E_rT)
pos = np.where(E_rT == minimo)
posf, = zip(pos[0], pos[1])
R_base = R_tot + (posf[0]-3)*dR
L_base = L_tot + (posf[1]-3)*dL

# =============================================================================
#     AJUSTE FINO DE R_tot y L_tot
# =============================================================================
i = 0
E_rT_fino = np.zeros((11, 11))
R_fino = np.zeros(11)
L_fino = np.zeros(11)
for R in np.arange(R_base-.5, R_base+.51, 0.1):
    j = 0
    for L in np.arange(L_base-5, L_base+5.1, 1):
        param, E_rT_i, T_c_def = star(M_tot, X, Y, R, L, T_c, num_capas)
        E_rT_fino[i, j] = E_rT_i
        L_fino[j] = L
        j += 1
    R_fino[i] = R
    i += 1

minimo_fino = np.min(E_rT_fino)
pos_fino = np.where(E_rT_fino == minimo_fino)
posf_fino, = zip(pos_fino[0], pos_fino[1])
R_final = R_fino[posf_fino[0]]
L_final = L_fino[posf_fino[1]]

# Volvemos a ejecutar el programa con los mejores valores encontrados.
param_f, E_rT_f, T_c_def = star(M_tot, X, Y, R_final, L_final, T_c, num_capas)

# =============================================================================
#     SALIDA DE DATOS
# =============================================================================
Capa = param_f['Capa']
Transporte = param_f['Transporte']
Ciclo = param_f['Ciclo']
Radio = param_f['Radio']
Presion = param_f['Presión']
Temperatura = param_f['Temperatura']
Masa = param_f['Masa']
Luminosidad = param_f['Luminosidad']
capas = len(Capa)

# -------------------------------------------------------------------------- #
# DATOS DE LA ESTRELLA
print('Parámetros de la estrella:')
print('    - Masa total: {:.1f}'.format(M_tot), 'x 10^33 g')
print('    - Radio total: {:.1f}'.format(round(R_final, 1)), 'x 10^10 cm')
print('    - Luminosidad total: {:.1f}'.format(int(L_final)), 'x 10^33 erg/s')
print('    - Temperatura central: {:.2f}'.format(T_c_def), 'x 10^7 K')
print('    - Fracción en masa de hidrógeno: {:.2f}'.format(X))
print('    - Fracción en masa de helio:: {:.2f}'.format(Y))
print('\n')

print('{:^100}'.format('Modelo numérico del interior estelar' '\n'))
print('{0:^5} | {1:^12} | {2:^7} | {3:^11} | {4:^11} | '
      '{5:^11} | {6:^11} | {7:^11}'.format(
        'Capa', 'Transporte', 'Ciclo', 'Radio', 'Presión',
        'Temperatura', 'Masa', 'Luminosidad'))
print('-'*100)
for k in range(capas):
    print('{0:^5.0f} | {1:^12} | {2:^7} | {3:>11.4f} | {4:>11.4f} | '
          '{5:>11.4f} | {6:>11.4f} | {7:>11.4f}'.format(
            Capa[k], Transporte[k], Ciclo[k], Radio[k], Presion[k],
            Temperatura[k], Masa[k], Luminosidad[k]))
print('\n')

# -------------------------------------------------------------------------- #
# MÍNIMO ERROR RELATIVO
print('Error relativo total mínimo obtenido: {:.2f}%'.format(
        round(E_rT_f, 4)*100), '\n')

fig2, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
fig2.suptitle('Mínimo error relativo total', fontsize=18)

X1, Y1 = np.mgrid[slice(R_init, R_fin + dR, dR),
                  slice(L_init, L_fin + dL, dL)]
Z1 = E_rT*100
Z1 = Z1[:-1, :-1]

X2, Y2 = np.mgrid[slice(R_base-.5, R_base+.6, 0.1),
                  slice(L_base-5, L_base+6, 1)]
Z2 = E_rT_fino*100
Z2 = Z2[:-1, :-1]

cmap = plt.get_cmap('viridis')

levels1 = MaxNLocator(nbins=100).tick_values(0, Z1.max())
norm1 = BoundaryNorm(levels1, ncolors=cmap.N, clip=True)

levels2 = MaxNLocator(nbins=100).tick_values(0, Z2.max())
norm2 = BoundaryNorm(levels2, ncolors=cmap.N, clip=True)

im1 = ax1.pcolormesh(X1, Y1, Z1, cmap=cmap, norm=norm1)
fig2.colorbar(im1, ax=ax1)
ax1.set_title('Ajuste de R_tot y L_tot')
ax1.set_xlabel('Radio (10^10 cm)')
ax1.set_ylabel('Luminosidad (10^33 erg/s)')

im2 = ax2.pcolormesh(X2, Y2, Z2, cmap=cmap, norm=norm2)
fig2.colorbar(im2, ax=ax2)
ax2.set_title('Ajuste fino')
ax2.set_xlabel('Radio (10^10 cm)')
ax2.set_ylabel('Luminosidad (10^33 erg/s)')

con = ConnectionPatch(xyA=(R_base + .05, L_base + .5),
                      xyB=(R_base + .25, L_base + 2.5),
                      coordsA="data", coordsB="data",
                      axesA=ax2, axesB=ax1, color='r',
                      arrowstyle="<|-", mutation_scale=20, linewidth=3)
ax2.add_artist(con)

fig2.tight_layout(rect=[0, 0.05, 1, 0.95])

plt.show()

print('\n')

# -------------------------------------------------------------------------- #
# PARÁMETROS EN FUNCIÓN DEL RADIO
fig, ax = plt.subplots(figsize=(10, 5))
ax.plot(Radio, Presion/max(Presion), 'r', label='Presión')
ax.plot(Radio, Temperatura/max(Temperatura), 'b', label='Temperatura')
ax.plot(Radio, Masa/max(Masa), 'g', label='Masa')
ax.plot(Radio, Luminosidad/max(Luminosidad), 'y', label='Luminosidad')
ax.set_title('Interior estelar', fontsize=18)
ax.set(xlabel='Radio (10^10 cm)', ylabel='Fracción del valor total (%)')
ax.legend()
