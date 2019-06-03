# -*- coding: utf-8 -*-


def ciclo(T_i, X, Y):
    """
    Ritmo de generación de energía

    Determinación del método y ritmo de generación de energía (si la hubiese)
    en una determinada capa de la estrella.

    Ciclos de generación de energía:
        PP: cadena protón-protón
        CN: ciclo CN
    Se calcula la generación de energía para ambos ciclos y se asume como
    único aquel con el que se obtiene mayor energía.
    En el caso de que la temperatura caiga fuera del intervalo asumido, la
    generación de energía será nula.

    Input:
        T_i: temperatura en la capa
        X: fracción en masa de hidrógeno en la estrella
        Y: fracción en masa de helio en la estrella

    Output:
        CICLO: ciclo de generación de energía (PP, CN o nulo)
        E1, nu: parámetros dependientes del ciclo de generación de energía y
                de la temperatura
    """

    # CADENA PROTÓN-PROTÓN
    if 0.4 < T_i <= 0.6:
        E1_pp = 10**(-6.84)
        nu_pp = 6.0
    elif 0.6 < T_i <= 0.95:
        E1_pp = 10**(-6.04)
        nu_pp = 5.0
    elif 0.95 < T_i <= 1.2:
        E1_pp = 10**(-5.56)
        nu_pp = 4.5
    elif 1.2 < T_i <= 1.65:
        E1_pp = 10**(-5.02)
        nu_pp = 4.0
    elif 1.65 < T_i <= 2.4:
        E1_pp = 10**(-4.40)
        nu_pp = 3.5
    else:
        E1_pp = 0.
        nu_pp = 0.
    # Ecuación para el cálculo del ritmo de generación de energía,
    # asumiendo cadena protón-protón.
    E_pp = E1_pp * (X**2) * (10*T_i)**nu_pp

    # CICLO CN
    if 1.2 < T_i <= 1.6:
        E1_CN = 10**(-22.2)
        nu_CN = 20.0
    elif 1.6 < T_i <= 2.25:
        E1_CN = 10**(-19.8)
        nu_CN = 18.0
    elif 2.25 < T_i <= 2.75:
        E1_CN = 10**(-17.1)
        nu_CN = 16.0
    elif 2.75 < T_i <= 3.60:
        E1_CN = 10**(-15.6)
        nu_CN = 15.0
    elif 3.60 < T_i <= 5.0:
        E1_CN = 10**(-12.5)
        nu_CN = 13.0
    else:
        E1_CN = 0.
        nu_CN = 0.
    # Ecuación para el cálculo del ritmo de generación de energía,
    # asumiendo ciclo CN.
    Z = 1 - X - Y  # Fracción en masa de elementos pesados.
    E_CN = (1/3) * E1_CN * X * Z * (10*T_i)**nu_CN

    # CICLO PREDOMINANTE
    # El ciclo que genere mayor energía será el predominante y, por tanto,
    # el que usaremos en nuestros cálculos.
    if E_pp > E_CN:
        CICLO = 'PP'
        E1 = E1_pp
        nu = nu_pp
    elif E_pp < E_CN:
        CICLO = 'CN'
        E1 = E1_CN
        nu = nu_CN
    else:
        CICLO = '--'
        E1 = 0.
        nu = 0.
    return(CICLO, E1, nu)
