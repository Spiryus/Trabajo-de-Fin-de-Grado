# -*- coding: utf-8 -*-


import numpy as np
import energia
import derydif as dd
import politropo as poli
import valor_estimado as vest
import valor_calculado as vcal
import valores_frontera as vf


def conveccion(T_c, rf, h, X, Y, k_prima):
    """
    Fase B del modelo numérico de interior estelar

    Integración de los distintos parámetros físicos del núcleo convectivo de la
    estrella (radio, presión, temperatura, masa y luminosidad), así como la
    generación de energía (si la hubiese) y la forma en que esta se transporta.
    Este cálculo se realiza para cada capa del núcleo estelar, empezando desde
    el centro de la estrella, con una estimación de la temperatura central,
    hasta la frontera radiativo-convectiva.

    Para iniciar la integración necesitamos conocer los valores de los
    parámetros en las tres primeras capas, para los cual utilizamos las
    mismas ecuaciones que en la Fase 0 del modelo completo ("estrella.py").

    Input:
        T_c: valor asumido para la temperatura central de la estrella
        rf: radio en la frontera ratiativo-convectiva
        h: paso de integración
            h = (radio inicial de integración) / (número de capas)
        X: fracción en masa de hidrógeno en la estrella
        Y: fracción en masa de helio en la estrella
        k_prima: constante politrópica, K'

    Output:
        param: narray
               parámetros para cada capa del núcleo de la estrella
        Frontera: narray
                  valores de los parámetros físicos en la frontera
        Ciclo: list
               ciclo de generación de energía en cada capa
    """

    Z = 1 - X - Y  # Fracción en masa de elementos pesados.
    mu = 1 / (2*X + (3/4)*Y + (1/2)*Z)  # Peso molecular medio.

    E_rel_max = 0.0001  # Error relativo máximo que aceptamos.

    # Creamos arrays para ir metiendo los valores de cada parámetro y sus
    # primeras derivadas.
    r = np.zeros(200)
    P = np.zeros(200)
    T = np.zeros(200)
    M = np.zeros(200)
    L = np.zeros(200)
    # dP = np.zeros(200)
    dT = np.zeros(200)
    dM = np.zeros(200)
    dL = np.zeros(200)
    Ciclo = []  # Ciclo de generación de energía, si lo hubiera.

# =============================================================================
#      CÁLCULO DE LAS TRES PRIMERAS CAPAS
# =============================================================================
    T_i = T_c  # Empezamos con la temperatura central asumida.

    for i in np.arange(3):
        # Ciclo de generación de energía.
        CICLO, E1, nu = energia.ciclo(T_i, X, Y)
        # Fracción en masa según ciclo predominante.
        if CICLO == 'PP':
            X1 = X
            X2 = X
        elif CICLO == 'CN':
            X1 = X
            X2 = (1/3)*(1-X-Y)
        else:
            X1 = 0.
            X2 = 0.

        # Magnitudes
        r_i = h * i
        T_i = T_c - (0.008207 * (mu**2) * k_prima * (T_c**1.5) * (r_i**2))
        P_i = k_prima * (T_i**2.5)
        L_i = (0.006150 * E1 * (X1*X2) * (10**nu) * (mu**2) * (k_prima**2)
               * (T_c**(3+nu)) * (r_i**3))
        M_i = 0.005077 * mu * k_prima * (T_c**1.5) * (r_i**3)
        r[i] = r_i
        P[i] = P_i
        T[i] = T_i
        M[i] = M_i
        L[i] = L_i
        Ciclo.append('--')

        # Primeras derivadas.
        if r_i != 0:
            f_i = dd.derivadas(r_i, P_i, T_i, M_i, L_i,
                               mu, X, Y, transfer='CONVECTIVE')
            # dP[i] = f_i[0]
            dT[i] = f_i[1]
            dM[i] = f_i[2]
            dL[i] = f_i[3]

# =============================================================================
#     FASE B.1.: Integración desde el centro. Núcleo convectivo.
#                Resto de capas.
# =============================================================================
    loop1 = True  # LOOP PARA LA CAPA
    while loop1:
        # Radio en la nueva capa.
        r_i1 = h*(i+1)
        # Temperatura estimada.
        T_est = vest.temperatura_estimada(
                T[i], dT[i], dT[i-1], dT[i-2], h)

        loop2 = True  # LOOP PARA LA TEMPERATURA
        while loop2:
            # Presión estimada.
            P_est = poli.presion(T_est, k_prima)
            # Masa calculada.
            M_cal, dM_i1 = vcal.masa_calculada(
                    r_i1, P_est, T_est, M[i], dM[i], dM[i-1], mu, h)
            # Temperatura calculada.
            T_cal, dT_i1 = vcal.temperatura_calculada(
                    r_i1, P_est, T[i], T_est, M_cal, L[i],
                    dT[i], dT[i-1], mu, h, X, Y,
                    transfer='CONVECTIVE')
            # ¿T_cal = T_est?
            E_rel = abs((T_cal-T_est)/T_cal)  # Error relativo.
            if E_rel < E_rel_max:
                loop2 = False
            else:
                T_est = T_cal

        # Presión calculada.
        P_cal = poli.presion(T_cal, k_prima)
        # Luminosidad calculada.
        L_cal, dL_i1, Ciclo_i1 = vcal.luminosidad_calculada(
                r_i1, P_cal, T_cal, L[i], dL[i], dL[i-1], mu, h, X, Y)

        # ¿Hemos llegado a la frontera radiativo-convectiva?
        if r_i1 > rf + h:
            loop1 = False
        else:
            # Damos por válidos los valores de esta capa.
            r[i+1] = r_i1
            P[i+1] = P_cal
            T[i+1] = T_cal
            M[i+1] = M_cal
            L[i+1] = L_cal
            # dP[i+1] = NULL
            dT[i+1] = dT_i1
            dM[i+1] = dM_i1
            dL[i+1] = dL_i1
            Ciclo.append(Ciclo_i1)
            # Pasamos a la capa siguiente.
            i += 1

    # Agrupamos los cinco parámetros en una sola variable.
    param = np.array([r, P, T, M, L])

# =============================================================================
#     VALORES EN LA FRONTERA
# =============================================================================
    Frontera = vf.frontera(rf, [r[i-1], r[i]],
                               [P[i-1], P[i]],
                               [T[i-1], T[i]],
                               [M[i-1], M[i]],
                               [L[i-1], L[i]])

# ----------------------------------------------------------------------------
    return(param, Frontera, Ciclo)
