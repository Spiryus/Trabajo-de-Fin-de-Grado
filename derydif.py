# -*- coding: utf-8 -*-


import primera_derivada as pd


def derivadas(r_i, P_i, T_i, M_i, L_i, mu, X, Y,
              P_cte=False, T_cte=False, M_cte=False, L_cte=False, transfer=[]):
    """
    Cálculo de las primeras derivadas para los cuatro parámetros físicos

    Para su cálculo utilizamos 'primera_derivada.py' que hace uso de
    las 4 ecuaciones fundamentales del interior estelar.

    Input:
        r_i: radio
        P_i: presión
        T_i: temperatura
        M_i: masa
        L_i: luminosidad
        + Si alguno se está considerando constante (ej. P_cte=True) ha de
          especificarse.
        + Por defecto se considera transporte radiativo. En caso de estar
          calculando los valores en un capa convectiva se ha de especificar
          (transfer = CONVECTIVE).

    Output:
        f: lista
           Primera derivada para los 4 parámetros físicos
    """

    dP = pd.presion(r_i, P_i, T_i, M_i, mu, P_cte)
    dT = pd.temperatura(r_i, P_i, T_i, M_i, L_i, mu, X, Y, T_cte, transfer)
    dM = pd.masa(r_i, P_i, T_i, mu, M_cte)
    dL, CICLO = pd.luminosidad(r_i, P_i, T_i, mu, X, Y, L_cte)
    f = [dP, dT, dM, dL]
    return(f)


def diferencias(f_i1, f_i, f_im1, h):
    """
    Método de diferencias

    Cálculo de diferencias a partir de la primera derivada de capas anteriores.
    Debido a la complejidad para calcular derivadas de orden superior,
    utilizamos el método de diferencias para el cálculo de siguientes capas.

    Input:
        f_i1: primera derivada en la capa que estamos calculando (i+1)
        f_i: primera derivada en la capa anterior (i)
        f_im1: primera derivada en la capa i-1
        h: paso de integración
            h = (radio inicial de integración) / (número de capas)

    Output:
        delta: lista
               delta1: diferencia de orden 1
               delta2: diferencia de orden 2
    """

    delta1 = h*f_i1 - h*f_i
    delta2 = h*f_i1 - 2*h*f_i + h*f_im1
    delta = [delta1, delta2]
    return(delta)
