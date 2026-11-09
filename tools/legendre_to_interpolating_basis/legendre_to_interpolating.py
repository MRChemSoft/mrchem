import sympy as sp
import numpy as np
import numpy.polynomial.legendre as leg
from numpy import sqrt


def Lj(order, x):
    """Returns the Legendre Basis function of a given order."""
    if (order == 0):
        return 1.0
    if (order == 1):
        return x
    return (2*(order-1)+1)/((order-1) + 1)*x*Lj(order-1, x) - \
           (order-1)/(order)*Lj(order-2, x)


def scaling_L(k, j, x):
    """Returns the Legendre Scaling function defined on x \in [0, 1]"""
    # k serves as a dummy variable made to give it the same signature as
    # scaling_I
    return sqrt(2.0*j + 1)*Lj(j, 2.0*x - 1)


def scaling_I(k, j, x):
    """Returns the Legendre Interpolating Scaling function defined
       on \in [0, 1]"""
    weight = leg.leggauss(k)[1]/2
    roots = leg.Legendre.basis(k, domain=[0, 1]).roots()
    tmp = 0
    for i in range(k):
        tmp = tmp+sqrt(weight[j])*scaling_L(k, i, roots[j])*scaling_L(k, i, x)
    return tmp


def make_transition_matrix(scaling, order):
    """Given the scaling function and the order this function
       generates the transition matrix Ps<-basis, where s is
       the standard basis.
    """
    x = sp.symbols('x')
    transition_matrix = np.zeros([order, order], dtype=np.float64)
    for j in range(order):
        try:
            pol = sp.poly(scaling(order, j, x))
            tmp = np.zeros(order, dtype=np.float64)
            pol_arr = np.flip(np.array(pol.all_coeffs(), dtype=np.float64), 0)
            for i in range(pol_arr.size):
                tmp[i] = pol_arr[i]
            transition_matrix[:, j] = tmp
        except:
            transition_matrix[0, 0] = 1.0
    return transition_matrix


def convert_to_interpolating(order, legendre_matrix, reverse=False):
    """Converts a Legendre operator matrix to a Legendre Interpolating operator
       matrix.

       Keyword arguments:
       reverse -- go from Legendre Interpolating to Legendre (default False)
    """
    P_s_i = make_transition_matrix(scaling_I, order)
    P_s_l = make_transition_matrix(scaling_L, order)

    P_i_s = np.linalg.inv(P_s_i)
    P_l_s = np.linalg.inv(P_s_l)

    P_i_l = np.dot(P_i_s, P_s_l)
    P_l_i = np.dot(P_l_s, P_s_i)

    if (not reverse):
        return np.dot(P_i_l, np.dot(legendre_matrix, P_l_i))
    else:
        return np.dot(P_l_i, np.dot(legendre_matrix, P_i_l))
