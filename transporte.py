# -*- coding: utf-8 -*-


def n_mas_1(P_i1, T_i1, dP_i1, dT_i1):
    """
    Cálculo del parámetro 'n+1'.

    Cálculo del parámetro 'n+1' (Novotny 1973, Capítulo 6, Sección 4.4.),
    utilizado para determinar cuándo la hipótesis del transporte radiativo
    deja de ser válida.
    El transporte por convección comenzará a ser importante a partir del
    momento en el que 'n+1 <= 2.5'.

    Input:
        P_i1: presión
        T_i1: temperatura
        dP_i1: primera derivada de la presión en la capa que se está calculando
        dT_i1: primera derivada de la temperatura en la capa

    Output:
        n_mas_1: parámetro n+1
    """

    n_mas_1 = (T_i1/P_i1) * (dP_i1/dT_i1)
    return(n_mas_1)
