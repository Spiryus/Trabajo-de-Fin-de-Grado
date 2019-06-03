# -*- coding: utf-8 -*-


import math
import numpy as np
import transporte
import derydif as dd
import politropo as poli
import valor_estimado as vest
import valor_calculado as vcal
import valores_frontera as vf
import faseB as fb
import matplotlib.pyplot as plt


def star(M_tot, X, Y, R_tot, L_tot, T_c, num_capas=100):
    """
    Modelo numérico de interior estelar

    Modelo numérico del interior estelar que nos divide la estrella en varias
    secciones (capas) y nos calcula para cada una de ellas los parámetros
    físicos de la estrella (radio, presión, temperatura, masa y luminosidad),
    así como la generación de energía (si la hubiese) y la forma en que esta
    se transporta.

    El programa consta de dos partes principales y un ajuste para encajarlas
    minimizando el error. Se tienen así las siguientes etapas de cálculo:
        - Fase 0: dado que, para evitar errores tenemos que empezar la
                  integración en un radio inferior, en esta fase calculamos
                  las capas que quedarían fuera, es decir, las que hay entre
                  el radio de inicio de integración y el radio total de la
                  estrella.
        - Fase A: esta es la fase de integración desde la superficie (envoltura
                  radiativa). Partimos de la masa total y la composición
                  química de la estrella e integramos para todos los parámetros
                  hasta llegar al centro de la estrella, pasando por las partes
                  de transporte radiativo y convectivo.
                  Además se subdivide esta fase en varias subfases tomando
                  alguno(s) de los parámetros como constante y contrastando
                  cuando esta aproximación deja de ser válida.
                  Para poder iniciar la integración necesitamos conocer los
                  valores de los parámetros en las tres primeras capas, para
                  los cual se utilizan las mismas ecuaciones que en la Fase 0.
        - Fase B: esta es la fase de integración desde el centro (núcleo
                  convectivo). Partimos de la temperatura central de la
                  estrella e integramos para todos los parámetros hasta llegar
                  a la frontera radiativo-convectiva (calculada en la fase A).
        - Ajuste de la temperatura central: dado que la temperatura central es
                  un dato que escogemos arbitrariamente, ocurre que en la
                  frontera radiativo-convectiva al empalmar las fases A y B,
                  tenemos incongruencias. Debido a esto, se realiza un ajuste
                  del valor de la temperatura central para conseguir que el
                  paso de una fase a otra sea fluido y con sentido.

    Para los ajustes se utiliza el error relativo total, que ha de ser mínimo:
        E_rel_total = sqrt((E_rel_presion)^2 + (E_rel_masa)^2 + ...)
            E_rel_param = (param_radiativ - param_convective) / param_radiativ

    Para la integración numérica utilizamos un método predictor-corrector, con
    el cual realizamos primero una estimación del valor del parámetro y
    después un bucle que nos permita ajustar dicho parámetro hasta que la
    diferencia entre estimación y corrección sea mínima.

    Input:
        M_tot: masa total de la estrella
        X: fracción en masa de hidrógeno en la estrella
        Y: fracción en masa de helio en la estrella
        R_tot: valor inicial del radio
        L_tot: valor inicial de la luminosidad
        T_c: valor inicial de la temperatura central
        num_capas: número de capas
                   por defecto: 100

    Output:
        agrupado: dictionary
                  diccionario con los datos de todas las capas para cada
                  parámetro de la estrella
        mínimo: mínimo error relativo total
        T_c_def: temperatura central tras ajuste
    """

    Z = 1 - X - Y  # Fracción en masa de elementos pesados.
    mu = 1 / (2*X + (3/4)*Y + (1/2)*Z)  # Peso molecular medio.

    # Para evitar errores de convergencia, vamos a empezar la integración
    # con un radio algo inferior al valor total.
    R_ini = 0.9 * R_tot

    h = -R_ini / num_capas  # Paso de integración.
    E_rel_max = 0.0001  # Error relativo máximo que aceptamos.

    # Constantes necesarias.
    A1 = 1.9022 * mu * M_tot
    A2 = 10.645 * math.sqrt(M_tot / (mu*Z*(1+X)*L_tot))

# =============================================================================
# =============================================================================
# #      FASE 0: Capas por debajo del radio inicial, R_ini.
# =============================================================================
# =============================================================================
    # Creamos arrays para ir metiendo los valores de cada parámetro.
    r_0 = np.zeros(50)
    P_0 = np.zeros(50)
    T_0 = np.zeros(50)
    M_0 = np.zeros(50)
    L_0 = np.zeros(50)
    Ciclo_0 = []  # Ciclo de generación de energía, si lo hubiera.

    i = 0
    loop1 = True
    while loop1:
        # Magnitudes
        r_i = R_ini - h*(i+1)
        T_i = A1 * (1/r_i - 1/R_tot)
        if T_i <= 0:
            P_i = 0
        else:
            P_i = A2 * T_i**4.25
        L_i = L_tot
        M_i = M_tot

        # ¿Hemos alcanzado el radio total?
        if r_i > R_tot:
            loop1 = False
        else:
            # Damos por válidos los valores de esta capa.
            r_0[i] = r_i
            P_0[i] = P_i
            T_0[i] = T_i
            M_0[i] = M_i
            L_0[i] = L_i
            Ciclo_0.append('--')
            # Pasamos a la capa siguiente.
            i += 1

# =============================================================================
# =============================================================================
# #      FASE A.: Integración desde la superficie.
# =============================================================================
# =============================================================================
    # Creamos arrays para ir metiendo los valores de cada parámetro y sus
    # primeras derivadas.
    r_A = np.zeros(200)
    P_A = np.zeros(200)
    T_A = np.zeros(200)
    M_A = np.zeros(200)
    L_A = np.zeros(200)
    dP_A = np.zeros(200)
    dT_A = np.zeros(200)
    dM_A = np.zeros(200)
    dL_A = np.zeros(200)
    Ciclo_A = []  # Ciclo de generación de energía, si lo hubiera.

# =============================================================================
#      CÁLCULO DE LAS TRES PRIMERAS CAPAS
# =============================================================================
    for i in range(3):
        # Magnitudes
        r_i = R_ini + h*i
        T_i = A1 * (1/r_i - 1/R_tot)
        P_i = A2 * T_i**4.25
        L_i = L_tot
        M_i = M_tot
        r_A[i] = r_i
        P_A[i] = P_i
        T_A[i] = T_i
        M_A[i] = M_i
        L_A[i] = L_i

        # Primeras derivadas.
        f_i = dd.derivadas(r_i, P_i, T_i, M_i, L_i, mu, X, Y,
                           M_cte='TRUE', L_cte='TRUE')
        dP_A[i] = f_i[0]
        dT_A[i] = f_i[1]
        dM_A[i] = f_i[2]
        dL_A[i] = f_i[3]

        Ciclo_A.append('--')

# =============================================================================
#      FASE A.1.1.: Integración desde la superficie. Envoltura radiativa.
#                   Masa y luminosidad constantes.
# =============================================================================
    loop1 = True  # LOOP PARA LA CAPA
    while loop1:
        # Radio en la nueva capa.
        r_i1 = R_ini + h*(i+1)
        # Presión y temperatura estimados.
        P_est = vest.presion_estimada(
                P_A[i], dP_A[i], dP_A[i-1], dP_A[i-2], h)
        T_est = vest.temperatura_estimada(
                T_A[i], dT_A[i], dT_A[i-1], dT_A[i-2], h)

        loop2 = True  # LOOP PARA LA TEMPERATURA
        while loop2:
            loop3 = True  # LOOP PARA LA PRESIÓN
            while loop3:
                # Presión calculada.
                P_cal, dP_i1 = vcal.presion_calculada(
                        r_i1, P_A[i], P_est, T_est, M_tot,
                        dP_A[i], dP_A[i-1], mu, h)
                # ¿P_cal = P_est?
                E_rel = abs((P_cal-P_est)/P_cal)  # Error relativo.
                if E_rel < E_rel_max:
                    loop3 = False
                else:
                    P_est = P_cal

            # Temperatura calculada.
            T_cal, dT_i1 = vcal.temperatura_calculada(
                    r_i1, P_cal, T_A[i], T_est, M_tot, L_tot,
                    dT_A[i], dT_A[i-1], mu, h, X, Y)
            # ¿T_cal = T_est?
            E_rel = abs((T_cal-T_est)/T_cal)  # Error relativo.
            if E_rel < E_rel_max:
                loop2 = False
            else:
                T_est = T_cal

        # Masa calculada.
        M_cal, dM_i1 = vcal.masa_calculada(
                r_i1, P_cal, T_cal, M_tot, dM_A[i], dM_A[i-1], mu, h)
        # ¿Sigue siendo válida la hipótesis de masa constante?
        E_rel = abs((M_cal-M_tot)/M_cal)  # Error relativo.
        if E_rel > E_rel_max:
            loop1 = False
        else:
            # Damos por válidos los valores de esta capa.
            r_A[i+1] = r_i1
            P_A[i+1] = P_cal
            T_A[i+1] = T_cal
            M_A[i+1] = M_tot
            L_A[i+1] = L_tot
            dP_A[i+1] = dP_i1
            dT_A[i+1] = dT_i1
            dM_A[i+1] = dM_i1
            # dL_A[i+1] = 0.
            Ciclo_A.append('--')
            # Pasamos a la capa siguiente.
            i += 1

# =============================================================================
#      FASE A.1.2.: Integración desde la superficie. Envoltura radiativa.
#                   Masa variable y luminosidad constante.
# =============================================================================
    loop1 = True  # LOOP PARA LA CAPA
    while loop1:
        # Radio en la nueva capa.
        r_i1 = R_ini + h*(i+1)
        # Presión y temperatura estimados.
        P_est = vest.presion_estimada(
                P_A[i], dP_A[i], dP_A[i-1], dP_A[i-2], h)
        T_est = vest.temperatura_estimada(
                T_A[i], dT_A[i], dT_A[i-1], dT_A[i-2], h)

        loop2 = True  # LOOP PARA LA TEMPERATURA
        while loop2:
            loop3 = True  # LOOP PARA LA PRESIÓN
            while loop3:
                # Masa calculada.
                M_cal, dM_i1 = vcal.masa_calculada(
                        r_i1, P_est, T_est, M_A[i],
                        dM_A[i], dM_A[i-1], mu, h)
                # Presión calculada.
                P_cal, dP_i1 = vcal.presion_calculada(
                        r_i1, P_A[i], P_est, T_est, M_cal,
                        dP_A[i], dP_A[i-1], mu, h)
                # ¿P_cal = P_est?
                E_rel = abs((P_cal-P_est)/P_cal)  # Error relativo.
                if E_rel < E_rel_max:
                    loop3 = False
                else:
                    P_est = P_cal

            # Temperatura calculada.
            T_cal, dT_i1 = vcal.temperatura_calculada(
                    r_i1, P_cal, T_A[i], T_est, M_cal, L_tot,
                    dT_A[i], dT_A[i-1], mu, h, X, Y)
            # ¿T_cal = T_est?
            E_rel = abs((T_cal-T_est)/T_cal)  # Error relativo.
            if E_rel < E_rel_max:
                loop2 = False
            else:
                T_est = T_cal

        # Luminosidad calculada.
        L_cal, dL_i1, Ciclo = vcal.luminosidad_calculada(
                r_i1, P_cal, T_cal, L_tot, dL_A[i], dL_A[i-1], mu, h, X, Y)
        # ¿Sigue siendo válida la hipótesis de luminosidad constante?
        E_rel = abs((L_cal-L_tot)/L_cal)  # Error relativo.
        if E_rel > E_rel_max:
            loop1 = False
        else:
            # Damos por válidos los valores de esta capa.
            r_A[i+1] = r_i1
            P_A[i+1] = P_cal
            T_A[i+1] = T_cal
            M_A[i+1] = M_cal
            L_A[i+1] = L_tot
            dP_A[i+1] = dP_i1
            dT_A[i+1] = dT_i1
            dM_A[i+1] = dM_i1
            dL_A[i+1] = dL_i1
            Ciclo_A.append('--')
            # Pasamos a la capa siguiente.
            i += 1

# =============================================================================
#      FASE A.1.3.: Integración desde la superficie. Envoltura radiativa.
#                   Masa y luminosidad variables.
# =============================================================================
    # Creamos un array para el parámetro n+1; nos hará falta más adelante.
    n_mas_1 = np.zeros(200)

    loop1 = True  # LOOP PARA LA CAPA
    while loop1:
        # Radio en la nueva capa.
        r_i1 = R_ini + h*(i+1)
        # Presión y temperatura estimados.
        P_est = vest.presion_estimada(
                P_A[i], dP_A[i], dP_A[i-1], dP_A[i-2], h)
        T_est = vest.temperatura_estimada(
                T_A[i], dT_A[i], dT_A[i-1], dT_A[i-2], h)

        loop2 = True  # LOOP PARA LA TEMPERATURA
        while loop2:
            loop3 = True  # LOOP PARA LA PRESIÓN
            while loop3:
                # Masa calculada.
                M_cal, dM_i1 = vcal.masa_calculada(
                        r_i1, P_est, T_est, M_A[i],
                        dM_A[i], dM_A[i-1], mu, h)
                # Presión calculada.
                P_cal, dP_i1 = vcal.presion_calculada(
                        r_i1, P_A[i], P_est, T_est, M_cal,
                        dP_A[i], dP_A[i-1], mu, h)
                # ¿P_cal = P_est?
                E_rel = abs((P_cal-P_est)/P_cal)  # Error relativo.
                if E_rel < E_rel_max:
                    loop3 = False
                else:
                    P_est = P_cal

            # Luminosidad calculada.
            L_cal, dL_i1, Ciclo_i1 = vcal.luminosidad_calculada(
                    r_i1, P_cal, T_est, L_A[i],
                    dL_A[i], dL_A[i-1], mu, h, X, Y)
            # Temperatura calculada.
            T_cal, dT_i1 = vcal.temperatura_calculada(
                    r_i1, P_cal, T_A[i], T_est, M_cal, L_cal,
                    dT_A[i], dT_A[i-1], mu, h, X, Y)
            # ¿T_cal = T_est?
            E_rel = abs((T_cal-T_est)/T_cal)  # Error relativo.
            if E_rel < E_rel_max:
                loop2 = False
            else:
                T_est = T_cal

        # Cálculo del parámetro 'n+1'.
        n_mas_1_i1 = transporte.n_mas_1(P_cal, T_cal, dP_i1, dT_i1)
        n_mas_1[i+1] = n_mas_1_i1
        # ¿Sigue siendo válida la hipótesis de transporte radiativo?
        if n_mas_1_i1 <= 2.5 or (n_mas_1_i1 > n_mas_1[i] and n_mas_1[i] != 0):
            loop1 = False
        else:
            # Damos por válidos los valores de esta capa.
            r_A[i+1] = r_i1
            P_A[i+1] = P_cal
            T_A[i+1] = T_cal
            M_A[i+1] = M_cal
            L_A[i+1] = L_cal
            dP_A[i+1] = dP_i1
            dT_A[i+1] = dT_i1
            dM_A[i+1] = dM_i1
            dL_A[i+1] = dL_i1
            Ciclo_A.append(Ciclo_i1)
            # Pasamos a la capa siguiente.
            i += 1

    i_f = i  # Capa de la frontera radiativo-convectiva.
    # Cálculo de la constante politrópica k'; invariante en todo el modelo.
    k_prima = poli.cte_k(P_cal, T_cal)

# =============================================================================
#     FASE A.2.: Integración desde la superficie. Núcleo convectivo.
# =============================================================================
    loop1 = True  # LOOP PARA LA CAPA
    while loop1:
        # Radio en la nueva capa.
        r_i1 = R_ini + h*(i+1)
        # Temperatura estimada.
        T_est = vest.temperatura_estimada(
                T_A[i], dT_A[i], dT_A[i-1], dT_A[i-2], h)

        loop2 = True  # LOOP PARA LA TEMPERATURA
        while loop2:
            # Presión estimada.
            P_est = poli.presion(T_est, k_prima)
            # Masa calculada.
            M_cal, dM_i1 = vcal.masa_calculada(
                    r_i1, P_est, T_est, M_A[i], dM_A[i], dM_A[i-1], mu, h)
            # Temperatura calculada.
            if r_i1 == 0:
                T_cal = T_est
            else:
                T_cal, dT_i1 = vcal.temperatura_calculada(
                        r_i1, P_est, T_A[i], T_est, M_cal, L_A[i],
                        dT_A[i], dT_A[i-1], mu, h, X, Y,
                        transfer='CONVECTIVE')
            # ¿T_cal = T_est?
            E_rel = abs((T_cal-T_est)/T_cal)  # Error relativo.
            if E_rel < E_rel_max:
                loop2 = False
            else:
                T_est = T_cal

        # Presión calculada.
        if r_i1 == 0:
            P_cal = P_est
        else:
            P_cal = poli.presion(T_cal, k_prima)
        # Luminosidad calculada.
        L_cal, dL_i1, Ciclo_i1 = vcal.luminosidad_calculada(
                r_i1, P_cal, T_cal, L_A[i],
                dL_A[i], dL_A[i-1], mu, h, X, Y)

        # ¿Hemos llegado al centro de la estrella?
        if r_i1 < 0:
            loop1 = False
        else:
            # Damos por válidos los valores de esta capa.
            r_A[i+1] = r_i1
            P_A[i+1] = P_cal
            T_A[i+1] = T_cal
            M_A[i+1] = M_cal
            L_A[i+1] = L_cal
            dP_A[i+1] = dP_i1
            dT_A[i+1] = dT_i1
            dM_A[i+1] = dM_i1
            dL_A[i+1] = dL_i1
            Ciclo_A.append(Ciclo_i1)
            # Pasamos a la capa siguiente.
            i += 1

# =============================================================================
#     VALORES EN LA FRONTERA
# =============================================================================
    # Valor del radio para el cual 'n+1 = 2.5'.
    rf = np.interp(2.5, [n_mas_1[i_f+1], n_mas_1[i_f]],
                        [r_A[i_f+1], r_A[i_f]])

    Frontera_A = vf.frontera(rf, [r_A[i_f+1], r_A[i_f]],
                                 [P_A[i_f+1], P_A[i_f]],
                                 [T_A[i_f+1], T_A[i_f]],
                                 [M_A[i_f+1], M_A[i_f]],
                                 [L_A[i_f+1], L_A[i_f]])


# =============================================================================
# =============================================================================
# #       FASE B.: Integración desde el centro.
# =============================================================================
# =============================================================================
    # Calculamos la fase B para una determinada temperatura central.
    param, Frontera_B, Ciclo_B = fb.conveccion(T_c, rf, -h, X, Y, k_prima)

# =============================================================================
#     AJUSTE DE LA TEMPERATURA CENTRAL, T_c.
# =============================================================================
    # Analizamos si hay que aumentar o disminuir nuestro valor de T_c.
    Increase_T_c = False
    Decrease_T_c = False
    if Frontera_B[-1] < Frontera_A[-1]:
        Increase_T_c = True
        incremento = 0.05
    elif Frontera_B[-1] > Frontera_A[-1]:
        Decrease_T_c = True
        incremento = -0.05

    loopT = True
    while loopT:
        param, Frontera_B, Ciclo_B = fb.conveccion(T_c, rf, -h, X, Y, k_prima)

        if ((Increase_T_c and Frontera_B[-1] > Frontera_A[-1]) or Frontera_B[-1] == 0) or \
           ((Decrease_T_c and Frontera_B[-1] < Frontera_A[-1]) or Frontera_B[-1] == 0):
            loopT = False
        else:
            T_c += incremento

    # Valores entre los que está comprendida la temperatura central buscada.
    T_c_l = T_c - incremento
    T_c_init = min([T_c_l, T_c])
    T_c_fin = max([T_c_l, T_c])

# =============================================================================
#     AJUSTE FINO DE T_c
# =============================================================================
    T_prueba = np.zeros(501)  # Lista para los sucesivos T_c.
    E_rT = np.zeros(501)  # Lista para los sucesivos errores relativos totales.

    j = 0
    for T_c in np.arange(T_c_init, T_c_fin, E_rel_max):
        param, Frontera_B, Ciclo_B = fb.conveccion(T_c, rf, -h, X, Y, k_prima)

        E_r = abs((Frontera_A-Frontera_B)/Frontera_A)
        E_rT_j = np.sqrt(E_r[1]**2 + E_r[2]**2 + E_r[3]**2 + E_r[4]**2)
        T_prueba[j] = T_c
        E_rT[j] = E_rT_j
        j += 1

    minimo = min(E_rT)
    pos, = np.where(E_rT == minimo)
    T_c_def = float(T_prueba[pos])

    # Parámetros finales para la fase B.
    param, Frontera_B, Ciclo_B = fb.conveccion(T_c_def, rf, -h, X, Y, k_prima)
    r_B = param[0]
    P_B = param[1]
    T_B = param[2]
    M_B = param[3]
    L_B = param[4]

# =============================================================================
# =============================================================================
# #       RESULTADO FINAL
# =============================================================================
# =============================================================================
    r_B = np.trim_zeros(r_B, trim='b')
    r_B = r_B[:-1]
    P_B = np.trim_zeros(P_B, trim='b')
    P_B = P_B[:-1]
    T_B = np.trim_zeros(T_B, trim='b')
    T_B = T_B[:-1]
    M_B = np.trim_zeros(M_B, trim='b')
    M_B = M_B[:-1]
    L_B = np.trim_zeros(L_B, trim='b')
    L_B = L_B[:-1]
    Ciclo_B = Ciclo_B[:-1]

    r_A = np.trim_zeros(r_A, trim='b')
    r_A = r_A[:i_f+1]
    P_A = np.trim_zeros(P_A, trim='b')
    P_A = P_A[:i_f+1]
    T_A = np.trim_zeros(T_A, trim='b')
    T_A = T_A[:i_f+1]
    M_A = np.trim_zeros(M_A, trim='b')
    M_A = M_A[:i_f+1]
    L_A = np.trim_zeros(L_A, trim='b')
    L_A = L_A[:i_f+1]
    Ciclo_A = Ciclo_A[:i_f+1]

    r_0 = np.trim_zeros(r_0, trim='b')
    P_0 = np.trim_zeros(P_0, trim='b')
    T_0 = np.trim_zeros(T_0, trim='b')
    M_0 = np.trim_zeros(M_0, trim='b')
    L_0 = np.trim_zeros(L_0, trim='b')
    # Ciclo_0 = Ciclo_0

    r = np.concatenate((r_B, np.flip(r_A), r_0), axis=None)
    r = np.round(r, 4)
    P = np.concatenate((P_B, np.flip(P_A), P_0), axis=None)
    P = np.round(P, 4)
    T = np.concatenate((T_B, np.flip(T_A), T_0), axis=None)
    T = np.round(T, 4)
    M = np.concatenate((M_B, np.flip(M_A), M_0), axis=None)
    M = np.round(M, 4)
    L = np.concatenate((L_B, np.flip(L_A), L_0), axis=None)
    L = np.round(L, 4)
    i = np.arange(len(r))  # Número total de capas.

    Ciclo = Ciclo_B + Ciclo_A[::-1] + Ciclo_0

    radiative = ['RADIATIVO'] * (len(r_0)+len(r_A))
    convective = ['CONVECTIVO'] * len(r_B)
    Transp = convective + radiative

    agrupado = {'Capa': i,
                'Transporte': Transp,
                'Ciclo': Ciclo,
                'Radio': r,
                'Presión': P,
                'Temperatura': T,
                'Masa': M,
                'Luminosidad': L}

# ------------- Grágica del ajuste de la temperatura central --------------- #
#    fig, ax = plt.subplots()
#    ax.plot(T_prueba, E_rT*100)
#    ax.plot(T_c_def, minimo*100, 'ro')
#    ax.set(title='Ajuste de temperatura central',
#           xlabel='T_c (10^7 K)', ylabel='Error relativo total (%)')
#
# ----------------------------------------------------------------------------
    return(agrupado, minimo, T_c_def)
