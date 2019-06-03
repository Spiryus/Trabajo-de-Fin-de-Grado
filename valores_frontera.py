# -*- coding: utf-8 -*-


import numpy as np


def frontera(rf, r, P, T, M, L):
    """
    Magnitudes en la prontera radiativo-convectiva.

    Intermpolación lineal para calcular los valores de la temperatura,
    presión, masa y luminosidad a un valor específico del radio ('rf').
    Los objetos 'r', 'P', 'T', 'M', 'L' han de contener al menos dos elementos
    (un valor en la zona radiatva y otro en la convectiva).

    El valor del radio en la frontera, rf, viene como input y output para una
    mejor visualización.

    Input:
        rf: radio en la frontera radiativo-convectiva
        r: array_like
           radio
        P: array_like
           presión
        T: array_like
           temperatura
        M: array_like
           masa
        L: array_like
           luminosidad

    Output:
        Frontera: narray
                  valores de los parámetros físicos en la frontera
    """

    Pf = np.interp(rf, r, P)
    Tf = np.interp(rf, r, T)
    Mf = np.interp(rf, r, M)
    Lf = np.interp(rf, r, L)
    Frontera = np.array([rf, Pf, Tf, Mf, Lf])
    return(Frontera)
