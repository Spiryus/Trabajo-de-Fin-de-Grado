# -*- coding: utf-8 -*-


import derydif as dd


def presion_estimada(P_i, dP_i, dP_im1, dP_im2, h):
    """
    Presión estimada

    Una estimación de la presión en la capa i+1, obtenida a partir de los
    valores de capas anteriores. Para su cálculo empleamos la primera derivada
    de la presión en capas anteriores y el método de diferencias.

    Input:
        P_i: presión en la capa inmediatamente anterior (capa i)
        dP_i: primera derivada de la presión en la capa i
        dP_im1: primera derivada de la presión en la capa i-1
        dP_im2: primera derivada de la presión en la capa i-2
        h: paso de integración
            h = (radio inicial de integración) / (número de capas)

    Output:
        P_est: presión estimada en la capa i+1
    """

    delta_i = dd.diferencias(dP_i, dP_im1, dP_im2, h)
    P_est = P_i + h*dP_i + (1/2)*delta_i[0] + (5/12)*delta_i[1]
    return(P_est)


def temperatura_estimada(T_i, dT_i, dT_im1, dT_im2, h):
    """
    Temperatura estimada

    Una estimación de la temperatura en la capa i+1, obtenida a partir de los
    valores de capas anteriores. Para su cálculo empleamos la primera derivada
    de la temperatura en capas anteriores y el método de diferencias.

    Input:
        T_i: temperatura en la capa inmediatamente anterior (capa i)
        dT_i: primera derivada de la temperatura en la capa i
        dT_im1: primera derivada de la temperatura en la capa i-1
        dT_im2: primera derivada de la temperatura en la capa i-2
        h: paso de integración
            h = (radio inicial de integración) / (número de capas)

    Output:
        T_est: temperatura estimada en la capa i+1
    """

    delta_i = dd.diferencias(dT_i, dT_im1, dT_im2, h)
    T_est = T_i + h*dT_i + (1/2)*delta_i[0] + (5/12)*delta_i[1]
    return(T_est)
