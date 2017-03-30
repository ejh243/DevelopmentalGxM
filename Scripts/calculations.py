'''
Perform calculations on matrices
'''
import numpy as np

def corr2_coeff(a, b):
    '''
    Row-wise correlation coefficient calculation for two 2D arrays
    http://stackoverflow.com/questions/30143417/computing-the-correlation-coefficient-between-two-multi-dimensional-arrays
    :param a: np.array shape N x T
    :param b: np.array shape M x T
    :return: np.array N x M in which each element is a correlation coefficient
    '''

    # Row wise mean of input arrays & subtract from input arrays themselves
    A_mA = a - a.mean(1)[:, None]
    B_mB = b - b.mean(1)[:, None]

    # Sum of squares across rows
    ssA = (A_mA**2).sum(1)
    ssB = (B_mB**2).sum(1)

    # Get the correlation coefficients
    return np.dot(A_mA, B_mB.T) / np.sqrt(np.dot(ssA[:, None], ssB[None]))
