import numpy as n
from numpy import random as r

def malmquist(n_stars, rho):
    """
    Purpose: To demonstrate malmquist bias

    Usage: malmquist(10000,200)

    Arguments:
            n_stars: the number of stars to distribute

            rho: the radius of observation shpere in [pc]
            returns mean absolute magnitude of uniformly distributed
            magnitude limited survey

    Written: Ryan A. Arneson, UCI, 5/2012
    """
    mean_M = 4.8
    sigma_M = .3
    m_limit = 10
    #uniformly distribute the stars in a sphere
    radius = rho*(r.uniform(0,1,n_stars))**(1.0/3.0)
    abs_M = r.normal(mean_M,sigma_M,n_stars)
    M_kept =list()
    narrow_mag = list()
    for i in range(0,n_stars):
        if (abs_M[i] + 5*n.log10(radius[i])-5) < m_limit:
            M_kept.append(abs_M[i])
            
    for i in range(0,n_stars):
        if (abs_M[i] + 5*n.log10(radius[i])-5) <= 10.0 and (abs_M[i] + 5*n.log10(radius[i])-5) >= 9.9:
            narrow_mag.append(radius[i])
    
    print "Mean M:", n.mean(abs_M)
    print "Observed mean M:", n.mean(M_kept)
    print "Mean M - observed mean M:",n.mean(M_kept)-n.mean(abs_M)

    print "True mean distance [pc]:", n.mean(narrow_mag)
    print "Assumed mean distance [pc]:", ((10**(10.2/5.0))+(10**(10.1/5.0)))*.5
    
    return 

    """
    EXAMPLE OF OUTPUT:

    >>>malm.malmquist(500000,200)
    
    Mean M: 4.80009418171
    Observed mean M: 4.67533844006
    Mean M - observed mean M: -0.124755741651
    True mean distance [pc]: 114.773675712
    Assumed mean distance [pc]: 107.18033721

    The deviation of the true mean of M from the observed mean M is
    ~ -1.38*.3**2 = -0.1242 and is accurate to the fouth decimal point.

    The true mean distance is ~7.6 pc greater than the assumed mean distance.
    """
