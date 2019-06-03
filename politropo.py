# -*- coding: utf-8 -*-


def cte_k(P_i1, T_i1):
    """
    Constante politrópica, K'.

    Cálculo de la constante politrópica, suponiendo transporte adiabático
    con un índice adiabático: gamma = 5/3.

    Input:
        P_i1: presión
        T_i1: temperatura

    Output:
        k_prima: constante politrópica
    """

    k_prima = P_i1 / (T_i1**2.5)
    return(k_prima)


def presion(T_i1, k):
    """Expresión politrópica para el cálculo de la presión

    Suponemos que el transporte convectivo se produce de forma adiabática:
        Polítropo: P = K' * T**(y/(y-1))
        Asumiendo índice adiabático, y = 5/3: P = K' * T**2.5

    Input:
        T_i1: temperatura
        k: constante politrópica

    Output:
        P_i1: presión
    """

    if T_i1 <= 0:
        P_i1 = 0
    else:
        P_i1 = k * (T_i1**2.5)
    return(P_i1)
