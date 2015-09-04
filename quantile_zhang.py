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
    y = y.flatten()
    if (np.isnan(np.min(y)))&(rmnans==True):
        y = y[np.logical_not(np.isnan(y))]
    elif (np.isnan(np.min(y)))&(rmnans==False)
        raise Exception('You must not have nans in percentile calculation')
    n = len(y)
    # j is the largest integer no greater than (p*(n+1))
    j = math.floor(p*(n+1))
    # f is the interpolation factor
    f = p*(n+1) - j
    # Sort the array
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
