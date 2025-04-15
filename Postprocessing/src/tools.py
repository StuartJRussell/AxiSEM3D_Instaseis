
"""
This file contains general helper functions that are not essential to be in the 
main class.
"""

# Import modules
import numpy as np
from scipy.special import erf

def calculate_STF(source, time, magnitude):
    """
    The slip rate and its derivative are required by Instaseis but are not 
    explicitly outputted by AxiSEM3D.
    
    For the GaussianSTF class we can calculate it by replicating the AxiSEM3D 
    source code. Relevent code, from the AxiSEM3D root directory, is in:
            /SOLVER/src/core/source/elem_src/stf/
    
    Note that this may not currently work for all source types.
    """
    
    # Assign source variables, self explanatory
    half_duration = source['source_time_function']['half_duration']
    decay_factor = source['source_time_function']['decay_factor']
    time_shift = source['source_time_function']['time_shift']
    derivative_order = ['ERF', 'GAUSSIAN', 'FIRST_DERIVATIVE', 'RICKER'
        ].index(source['source_time_function']['use_derivative_integral']) - 1
    
    # STF
    STF = Gaussian_STF(time, half_duration, decay_factor, 
                                        time_shift, derivative_order)
    
    # Derivative of STF
    dSTF = Gaussian_STF(time, half_duration, decay_factor, 
                                        time_shift, derivative_order + 1)
    
    # STF should be normalised such that the absolute area is the source magnitude
    area_STF = np.trapz(np.abs(STF), x=time)
    STF = magnitude * STF / area_STF
    
    # Area of dSTF should be the sum of STF
    area_dSTF = np.trapz(np.abs(dSTF), x=time)
    dSTF = np.sum(np.abs(STF)) * dSTF / area_dSTF
    
    return STF, dSTF



def Gaussian_STF(time, half_duration, decay_factor, time_shift,
                                                 derivative_order):
    """
    This function simply calculates the Gaussian as required for the STF.
    Based on the code from AxiSEM3D.
    """

    # Some variables
    ispi = 1. / np.sqrt(np.pi)
    f = decay_factor / half_duration
    t = time - time_shift
    ft = f * t
    
    # Calculate curve based on parameters
    # ERF
    if derivative_order == -1:
        value = erf(ft) * 0.5 + 0.5
    # Gaussian
    elif derivative_order == 0:
        value = ispi * f * np.exp(-(ft ** 2.))
    # First derivative of Gaussian
    elif derivative_order == 1:
        value = ispi * (f ** 2.) * np.exp(-(ft ** 2.)) * (-2. * ft)
    # Second derivative of Gaussian
    elif derivative_order == 2:
        value = ispi * (f ** 3.) * np.exp(-(ft ** 2.)) * (4. * ft * ft - 2.)
    # Third derivative of Gaussian
    elif derivative_order == 3:
        value = (ispi * (f ** 4.) * np.exp(-(ft ** 2.)) * 
                            (12. * ft - 8. * ft ** 3.))
    
    return value



def calculate_GLL(n):
    """
    Calculates the GLL points.
    Finds the roots of the first derivative of a Legendre polynomial
    with order n-1.
    """
    
    # Order must be at least 2
    if n < 2:
        raise ValueError('N must be greater than 1')
    
    if n >= 2:
    
        # Define a large array between -1 and 1
        x = np.linspace(-1, 1, 100001)[1:-1]
        
        # Derivative of Legendre with n-1
        dP = dLegendre(n - 1, x)
        
        # Where this cross the x-axis are the GLL points
        indices = np.where((np.sign(dP[1:] * dP[:-1])==-1.))
        
        # GLL points
        gll = np.concatenate((np.array([-1.]), np.sort(x[indices]), np.array([1.])))
        
        # If odd then insert zero
        if n % 2 != 0:
            gll = np.insert(gll, int(len(gll)/2), 0)

    return gll



def calculate_GLJ(n):
    """
    Calculates the GLJ points.
    Finds the roots of the first derivative of a Jacobi polynomial
    with order n.
    """
    
    # Define a large array between -1 and 1
    x = np.linspace(-1, 1, 100001)[1:-1]
    
    # Get Jacobi
    y = Jacobi(n, x)
    
    # Differentiate it
    dJ = np.diff(y) / np.diff(x)
    
    # Where this cross the x-axis are the GLL points
    indices = np.where((np.sign(dJ[1:] * dJ[:-1]) == -1.))
    
    # GLL points
    glj = np.concatenate((np.array([-1.]), np.sort(x[indices]), np.array([1.])))

    return glj



def Legendre(n, x):
    """
    Calculate Legendre polynomials of order n at points x.
    Use Bonnet's recursion formula.
    """
    
    # We are required to explicitly define n=0 and n=1
    # They are easy at least
    if n == 0:
        y = np.ones(x.shape)
    elif n == 1:
        y = x
    
    # For n>1 we use recursion
    elif n >= 2:
        Pn1 = Legendre(n - 1, x)
        Pn2 = Legendre(n - 2, x)
        y = ((2. * (n - 1.) + 1.) * x * Pn1 - (n - 1.) * Pn2) / n
    
    return y



def dLegendre(n, x):
    """
    Calculate first derivative of Legendre polynomial of order n at points x.
    """
    
    # Exceptions for lower order polynomials
    if n == 0:
        y = np.zeros(x.shape)
    
    elif n >= 1:
        Pn = Legendre(n, x)
        Pn1 = Legendre(n - 1, x)
        y = n * (Pn1 - x * Pn) / (1. - x**2.)
    
    return y



def Jacobi(n, x):
    """
    Calculate Jacobi polynomials of order n at points x.
    This is equation A27 of Nissen-Meyer et al. 2007.
    """
    
    if n == 0:
        raise ValueError('N must be greater than 0')
    
    else:
        
        # Get Legendre polynomials
        Pn = Legendre(n, x)
        Pn1 = Legendre(n + 1, x)
        y = (Pn + Pn1) / (1. + x)
    
    return y



def calculate_G1(n, x):
    """
    The spatial derivatives of the basis functions, valid for eta direction of axial
    elements.
    
    This is equation A27 of Nissen-Meyer et al. 2007. Note that there is a typo in
    the seventh part of this equation in the paper, and the `I' in the demonenator
    should be `i'.
    """
    
    G = np.zeros((n+1, n+1))
    for i in range(n+1):
        for j in range(n+1):
            
            if i == 0 and j == 0:
                G[j, i] = -(n * (n + 2.)) / 6.
                
            elif i == 0 and 1 <= j <= (n - 1):
                G[j, i] = ((2. * ((-1.) ** n) * Jacobi(n, x[j])) /
                                ((1. + x[j]) * (n + 1.)))
                
            elif i == 0 and j == n:
                G[j, i] = ((-1.) ** n) / (n + 1.)
                
            elif 1 <= i <= (n - 1) and j == 0:
                G[j, i] = ((((-1.) ** (n + 1.)) * (n + 1.)) /
                            (2. * Jacobi(n, x[i]) * (1. + x[i])))
                
            elif 1 <= i <= (n - 1) and 1 <= j <= (n - 1) and i != j:
                G[j, i] = (1. / (x[j] - x[i])) * (Jacobi(n, x[j]) /
                                                  Jacobi(n, x[i]))
                
            elif 1 <= i <= (n - 1) and i == j:
                G[j, i] = -1. / (2. * (1. + x[i]))
            
            elif 1 <= i <= (n - 1) and j == n:
                G[j, i] = 1. / ((Jacobi(n, x[i]) * (1. - x[i])))
                
            elif i == n and j == 0:
                G[j, i] = (((-1.) ** (n + 1.)) * (n + 1.)) / 4.
                
            elif i == n and 1 <= j <= (n - 1):
                G[j, i] = (-1.) * Jacobi(n, x[j]) / (1. - x[j])
            
            elif i == n and j == n:
                G[j, i] = (n * (n + 2.) - 1.) / 4.
    
    return G



def calculate_G2(n, gll):
    """
    The spatial derivatives of the basis functions, valid for both directions of
    non-axial elements and for eta direction of axial ones.
    
    This is equation A20 of Nissen-Meyer et al. 2007.
    """
    
    G = np.zeros((n + 1, n + 1))
    for i in range(n + 1):
        for j in range(n + 1):
        
            if i == 0 and j == 0:
                G[j, i] = -(n * (n + 1.)) / 4.
    
            elif i == n and j == n:
                G[j, i] = (n * (n + 1.)) / 4.

            elif i == j:
                G[j, i] = 0.
            
            else:
                G[j, i] = ((Legendre(n, gll[j]) / Legendre(n, gll[i])) * 
                                (1. / (gll[j] - gll[i])))

    return G
