def quantile_R(x, p, itype=7, fraction=False, rmnans=False):
    """quantile function used in R

    Calculates quantiles in the same way as the quantile function in R.
    There are nine interpolation methods. type 1-3 are discrete methods and
    4-9 are continuous methods. 4 is a basic linear interpolation. 5 is the
    default method used in Matlab. 7 is the default used in R and Python.
    See: Hyndman, R. J. and Fan, Y. (1996) Sample quantiles in statistical
    packages, American Statistician, 50, 361-365.
    
    Arguments
    x -- array
    p -- percentile
    itype -- interpolation method
    fraction -- boolean indicates if percentile is a fraction 0<p<1.
    rmnans -- boolean indicates whether or not to remove nans

    Returns
    q -- quantile at pth percentile
    """
    import numpy as np
    import math
    # Convert the percentile to a fraction
    if not fraction:
        p = p/100.
    # Flatten the array
    x = x.flatten()
    # remove nans if desired
    if (np.isnan(np.min(x)))&(rmnans==True):
        x = x[np.logical_not(np.isnan(x))]
    elif (np.isnan(np.min(x)))&(rmnans==False):
        raise Exception('You must not have nans in percentile calculation')
    # Get number of samples
    n = len(x)
    # Sort
    x = np.sort(x)
    # Switch case functions for interpolation type
    def one(x,p,n):
        m = 0.
        j = np.floor(p*n + m)
        g = n*p + m - j
        if g==0: 
            gamma = 0.
        else:
            gamma = 1.
        return (1. - gamma)*x[j-1] + gamma*x[j]
    
    def two(x,p,n):
        m = 0.
        j = np.floor(p*n + m)
        g = n*p + m - j
        if g==0:
            gamma = 0.5
        else:
            gamma = 1.
        return (1. - gamma)*x[j-1] + gamma*x[j]

    def three(x,p,n):
        m = -0.5
        j = np.floor(p*n + m)
        g = n*p + m - j
        if (g==0)&(j%2==0):
            gamma = 0.
        else:
            gamma = 1.
        return (1. - gamma)*x[j-1] + gamma*x[j]

    def four(x,p,n):
        m = 0.
        j = np.floor(p*n + m)
        gamma = n*p + m - j
        return (1. - gamma)*x[j-1] + gamma*x[j]
    
    def five(x,p,n):
        m = 0.5
        j = np.floor(p*n + m)
        gamma = n*p + m - j
        return (1. - gamma)*x[j-1] + gamma*x[j]
    
    def six(x,p,n):
        m = p
        j = np.floor(p*n + m)
        gamma = n*p + m - j
        return (1. - gamma)*x[j-1] + gamma*x[j]
        
    def seven(x,p,n):
        m = 1. - p
        j = np.floor(p*n + m)
        gamma = n*p + m - j
        return (1. - gamma)*x[j-1] + gamma*x[j]

    def eight(x,p,n):
        m = (p+1)/3.
        j = np.floor(p*n + m)
        gamma = n*p + m - j
        return (1. - gamma)*x[j-1] + gamma*x[j]

    def nine(x,p,n):
        m = p/4. + 3./8.
        j = np.floor(p*n + m)
        gamma = n*p + m - j
        return (1. - gamma)*x[j-1] + gamma*x[j]

    switcher = {1: one, 2: two, 3: three, 4: four, 5: five, 6: six,
            7: seven, 8: eight, 9: nine}
    return switcher[itype](x,p,n)

def quantile_zhang(y, p, fraction=False, rmnans=False):
    """Caclulate the pth percentile value of an array y using the Zhang method.
    
    The linear interpolation method used is outlined in Zhang et al., 2005,
    Avoiding Inhomogeneity in Percentile-Based Indices of Temperature Extremes,
    Journal of Climate, vol. 18. The interpolation is given as:
        Qp =  (1-f)*y_j + f*y_j+1
    p<1/(n+1) are set to the smallest value in y
    p>n/(n+1) are set to the largest value in y
    
    Arguments
    y -- array
    p -- pth percentile
    fraction -- boolean indicates if percentile is a fraction 0<p<1.
    rmnans -- boolean indicates whether or not to remove nans
    
    Returns
    Qp -- qualtile of pth percentile
    """
    import numpy as np
    import math
    # Convert the percentile to a fraction
    if not fraction:
        p = p/100.
    # Flatten array
    #y = y.flatten()
    if (np.isnan(np.min(y)))&(rmnans==True):
        y = y[np.logical_not(np.isnan(y))]
    elif (np.isnan(np.min(y)))&(rmnans==False):
        raise Exception('You must not have nans in percentile calculation')
    n = y.shape[0]
    # j is the largest integer no greater than (p*(n+1))
    j = math.floor(p*(n+1))
    # f is the interpolation factor
    f = p*(n+1) - j
    if y.ndim==3:
        Qp = np.ones((y.shape[1],y.shape[2]))*np.nan
        for lat in range(y.shape[1]):
            for lon in range(y.shape[2]):
                # Sort the array
                array = np.sort(y[:,lat,lon])
                if j>=n:
                    Qp[lat,lon] = array[-1]
                elif j<1:
                    Qp[lat,lon] = array[0]
                else:
                    yj = array[j-1]
                    yjp = array[j]
                    Qp[lat,lon] = (1-f)*yj + f*yjp
    else:
        y = np.sort(y)
        if j>=n:
            Qp = y[-1]
        elif j<1:
            Qp = y[0]
        else:
            yj = y[j-1]
            yjp = y[j]
            Qp = (1-f)*yj + f*yjp
    return Qp

def quantile_climpact(y,p,fraction=False):
    import numpy as np
    import math
    import sys
    if not fraction:
        p = p/100.
    n = y.shape[0]
    a, b = 1.0/3.0, 1.0/3.0
    nppm = a + p*(n + 1 - a - b) - 1
    fuzz = 4*sys.float_info.epsilon
    j = math.floor(nppm + fuzz)
    if (abs(nppm - j)<=fuzz):
        h = 0
    else:
        h = nppm - j
    right_elem = max(0, min(int(j) + 1, n - 1))
    left_elem = max(0, min(int(j), n - 1))
    if y.ndim==3:
        Qp = np.ones((y.shape[1],y.shape[2]))*np.nan
        for lat in range(y.shape[1]):
            for lon in range(y.shape[2]):
                # Sort the array
                array = np.sort(y[:,lat,lon])
                if h==1:
                    Qp[lat,lon] = array[right_elem]
                elif h==0:
                    Qp[lat,lon] = array[left_elem]
                else:
                    Qp[lat,lon] = ((1 - h)*array[left_elem] + h*array[right_elem])
    else:
        y = np.sort(y)
        if h==1:
            Qp = y[right_elem]
        elif h==0:
            Qp = y[left_elem]
        else:
            Qp = (1 - h)*y[left_elem] + h*y[right_elem]
    return Qp
