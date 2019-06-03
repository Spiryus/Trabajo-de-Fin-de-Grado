# -*- coding: utf-8 -*-


"""
Cálculo de la primera derivada con respecto del radio

Cálculo de la primera derivada con respecto del radio para los distintos
parámetros físicos: presión, temperatura, masa y luminosidad.
 - Si alguno de los parámetros se está considerando constante (ej. P_cte=True),
   su derivada será nula.
 - Exceptuando la temperatura, para el resto de parámetros las ecuaciones
   en los casos radiativo y convectivo coinciden.
"""

import energia


def presion(r_i, P_i, T_i, M_i, mu, P_cte=False):
    """
    Primera derivada de la presión

    Input:
        r_i: radio
        P_i: presión
        T_i: temperatura
        M_i: masa
        mu: peso molecular medio
            mu = 1 / (2X + (3/4)Y + (1/2)Z)

    Output:
        dP_i: primera derivada de la presión en r_i
    """

    if P_cte:
        dP_i = .0
    else:
        # Ecuación de equilibrio hidrostático.
        dP_i = -8.084 * mu * (P_i/T_i) * (M_i/(r_i**2))
    return(dP_i)


def temperatura(r_i, P_i, T_i, M_i, L_i, mu, X, Y, T_cte=False, transfer=[]):
    """
    Primera derivada de la temperatura

    Input:
        r_i: radio
        P_i: presión
        T_i: temperatura
        M_i: masa
        L_i: luminosidad
        mu: peso molecular medio
            mu = 1 / (2X + (3/4)Y + (1/2)Z)
        X: fracción en masa de hidrógeno en la estrella
        Y: fracción en masa de helio en la estrella
        + Por defecto se considera transporte radiativo. En caso de estar
          calculando los valores en un capa convectiva se ha de especificar
          (transfer = CONVECTIVE).

    Output:
        dT_i: primera derivada de la temperatura en r_i
    """

    if T_cte:
        dT_i = .0
    elif transfer == 'CONVECTIVE':
        # Ecuación de transporte de energía (caso convectivo).
        dT_i = -3.234 * mu * (M_i/(r_i**2))
    else:
        Z = 1 - X - Y  # Fracción en masa de elementos pesados.
        # Ecuación de transporte de energía (caso radiativo).
        dT_i = (-0.01679 * Z * (X+1) * (mu**2)
                * ((P_i**2)/(T_i**8.5)) * (L_i/(r_i**2)))
    return(dT_i)


def masa(r_i, P_i, T_i, mu, M_cte=False):
    """
    Primera derivada de la masa

    Input:
        r_i: radio
        P_i: presión
        T_i: temperatura
        mu: peso molecular medio
            mu = 1 / (2X + (3/4)Y + (1/2)Z)

    Output:
        dM_i: primera derivada de la masa en r_i
    """

    if M_cte:
        dM_i = .0
    else:
        # Ecuación de continuidad de la masa.
        dM_i = 0.01523 * mu * (P_i/T_i) * r_i**2
    return(dM_i)


def luminosidad(r_i, P_i, T_i, mu, X, Y, L_cte=False):
    """
    Primera derivada de la luminosidad

    Input:
        r_i: radio
        P_i: presión
        T_i: temperatura
        mu: peso molecular medio
            mu = 1 / (2X + (3/4)Y + (1/2)Z)
        X: fracción en masa de hidrógeno en la estrella
        Y: fracción en masa de helio en la estrella

    Output:
        dL_i: primera derivada de la luminosidad en r_i
        CICLO: ciclo de generación de energía (PP, CN o nulo)
    """

    if L_cte:
        dL_i = .0
        CICLO = '--'
    else:
        # Ritmo de generación de energía.
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
        # Ecuación de equilibrio energético.
        dL_i = (0.01845 * E1 * (X1*X2) * (10**nu) * (mu**2)
                * (P_i**2) * (T_i**(nu-2)) * (r_i**2))
    return(dL_i, CICLO)
