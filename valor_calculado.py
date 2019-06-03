# -*- coding: utf-8 -*-


"""
Valor calculado para los 4 parámetros físicos

Valor calculado en la nueva capa para los 4 parámetros físicos: presión,
temperatura, masa y luminosidad.
Para su cálculo se emplean valores calculados o estimados en la nueva capa
para alguno(s) de los parámetros fisicos, lo que nos permite obtener la
primera derivada en la nueva capa del parámetro que estamos calculando, y
el método de diferencias.
 - Si alguno de los parámetros se está considerando constante (ej. P_cte=True),
   su derivada será nula.
 - Exceptuando la temperatura, para el resto de parámetros las ecuaciones
   en los casos radiativo y convectivo coinciden.
"""

import primera_derivada as pd
import derydif as dd


def presion_calculada(r_i1, P_i, P_i1, T_i1, M_i1, dP_i, dP_im1,
                      mu, h, P_cte=False):
    """
    Presión calculada

    Input:
        r_i1: radio en la nueva capa (la que estamos calculando, i+1)
        P_i: presión en la capa inmediatamente anterior (capa i)
        P_i1: una estimación de la presión para la nueva capa
        T_i1: una estimación de la temperatura en la nueva capa
        M_i1: una estimación de la masa en la nueva capa
        dP_i: primera derivada de la presión en la capa anterior (capa i)
        dP_im1: primera derivada de la presión en la capa i-1
        mu: peso molecular medio
            mu = 1 / (2X + (3/4)Y + (1/2)Z)
        h: paso de integración
            h = (radio inicial de integración) / (número de capas)

    Output:
        P_cal: presión calculada para la capa i+1
        dP_i1: primera derivada de la presión en la capa i+1
    """

    # Primera derivada en la nueva capa 'i+1'.
    dP_i1 = pd.presion(r_i1, P_i1, T_i1, M_i1, mu, P_cte)
    # Primeras diferencias en la nueva capa 'i+1'.
    delta_i1 = dd.diferencias(dP_i1, dP_i, dP_im1, h)

    P_cal = P_i + h*dP_i1 - (1/2)*delta_i1[0]
    return(P_cal, dP_i1)


def temperatura_calculada(r_i1, P_i1, T_i, T_i1, M_i1, L_i1, dT_i, dT_im1,
                          mu, h, X, Y, T_cte=False, transfer=[]):
    """
    Temperatura calculada

    Input:
        r_i1: radio en la nueva capa (la que estamos calculando, i+1)
        P_i1: una estimación de la presión en la nueva capa
        T_i: temperatura en la capa inmediatamente anterior (capa i)
        T_i1: una estimación de la temperatura en la nueva capa
        M_i1: una estimación de la masa en la nueva capa
        L_i1: una estimación de la luminosidad en la nueva capa
        dT_i: primera derivada de la temperatura en la capa anterior (capa i)
        dT_im1: primera derivada de la temperatura en la capa i-1
        mu: peso molecular medio
            mu = 1 / (2X + (3/4)Y + (1/2)Z)
        h: paso de integración
            h = (radio inicial de integración) / (número de capas)
        X: fracción en masa de hidrógeno en la estrella
        Y: fracción en masa de helio en la estrella
        + Por defecto se considera transporte radiativo. En caso de estar
          calculando los valores en un capa convectiva se ha de especificar
          (transfer = CONVECTIVE).

    Output:
        T_cal: temperatura calculada para la capa i+1
        dT_i1: primera derivada de la temperatura en la capa i+1
    """

    # Primera derivada en la nueva capa 'i+1'.
    dT_i1 = pd.temperatura(r_i1, P_i1, T_i1, M_i1, L_i1,
                           mu, X, Y, T_cte, transfer)
    # Primeras diferencias en la nueva capa 'i+1'.
    delta_i1 = dd.diferencias(dT_i1, dT_i, dT_im1, h)

    T_cal = T_i + h*dT_i1 - (1/2)*delta_i1[0]
    return(T_cal, dT_i1)


def masa_calculada(r_i1, P_i1, T_i1, M_i, dM_i, dM_im1,
                   mu, h, M_cte=False):
    """
    Masa calculada

    Input:
        r_i1: radio en la nueva capa (la que estamos calculando, i+1)
        P_i1: una estimación de la presión para la nueva capa
        T_i1: una estimación de la temperatura en la nueva capa
        M_i: masa en la capa inmediatamente anterior (capa i)
        dM_i: primera derivada de la masa en la capa anterior (capa i)
        dM_im1: primera derivada de la masa en la capa i-1
        mu: peso molecular medio
            mu = 1 / (2X + (3/4)Y + (1/2)Z)
        h: paso de integración
            h = (radio inicial de integración) / (número de capas)

    Output:
        M_cal: masa calculada para la capa i+1
        dM_i1: primera derivada de la masa en la capa i+1
    """

    # Primera derivada en la nueva capa 'i+1'.
    dM_i1 = pd.masa(r_i1, P_i1, T_i1, mu, M_cte)
    # Primeras diferencias en la nueva capa 'i+1'.
    delta_i1 = dd.diferencias(dM_i1, dM_i, dM_im1, h)

    M_cal = M_i + h*dM_i1 - (1/2)*delta_i1[0]
    return(M_cal, dM_i1)


def luminosidad_calculada(r_i1, P_i1, T_i1, L_i, dL_i, dL_im1,
                          mu, h, X, Y, L_cte=False):
    """
    Luminosidad calculada

    Input:
        r_i1: radio en la nueva capa (la que estamos calculando, i+1)
        P_i1: una estimación de la presión en la nueva capa
        T_i1: una estimación de la temperatura en la nueva capa
        M_i: luminosidad en la capa inmediatamente anterior (capa i)
        dL_i: primera derivada de la luminosidad en la capa anterior (capa i)
        dL_im1: primera derivada de la luminosidad en la capa i-1
        mu: peso molecular medio
            mu = 1 / (2X + (3/4)Y + (1/2)Z)
        h: paso de integración
            h = (radio inicial de integración) / (número de capas)
        X: fracción en masa de hidrógeno en la estrella
        Y: fracción en masa de helio en la estrella

    Output:
        L_cal: luminosidad calculada para la capa i+1
        dL_i1: primera derivada de la luminosidad en la capa i+1
        CICLO: ciclo de generación de energía (PP, CN o nulo)
    """

    # Primera derivada y ciclo de generación de energía predominante
    # en la nueva capa 'i+1'.
    dL_i1, CICLO = pd.luminosidad(r_i1, P_i1, T_i1, mu, X, Y, L_cte)
    # Primeras diferencias en la nueva capa 'i+1'.
    delta_i1 = dd.diferencias(dL_i1, dL_i, dL_im1, h)

    L_cal = L_i + h*dL_i1 - (1/2)*delta_i1[0] - (1/12)*delta_i1[1]
    return(L_cal, dL_i1, CICLO)
