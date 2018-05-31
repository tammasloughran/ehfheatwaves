#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
qtiler.py module contains functions that calculate quantiles using various
interpolation methods.

@author: Tammas Loughran
"""
import numpy as np
import math
import sys

class InvalidPercentileError(Exception):
    """Exception to be raised if the provided percentile is invalid."""
    def __init__(self, p):
        print(p, " is not a valid fraction value.")


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

    # Convert the percentile to a fraction
    if not fraction:
        p = p/100.
    if (p>1) or (p<0): raise InvalidPercentileError(p)
    # Flatten the array
    x = x.flatten()
    # remove nans if desired
    if (np.isnan(np.min(x)))&(rmnans==True):
        x = x[np.logical_not(np.isnan(x))]
    elif (np.isnan(np.min(x)))&(rmnans==False):
        raise Exception('You must not have nans in percentile calculation')
    if (p==1): return max(x)
    # Get number of samples
    n = len(x)
    # Sort
    x = np.sort(x)
    # Switch case functions for interpolation type
    def one(x,p,n):
        m = 0.
        j = int(np.floor(p*n + m))
        g = n*p + m - j
        if g==0:
            gamma = 0.
        else:
            gamma = 1.
        return (1. - gamma)*x[j-1] + gamma*x[j]

    def two(x,p,n):
        m = 0.
        j = int(np.floor(p*n + m))
        g = n*p + m - j
        if g==0:
            gamma = 0.5
        else:
            gamma = 1.
        return (1. - gamma)*x[j-1] + gamma*x[j]

    def three(x,p,n):
        m = -0.5
        j = int(np.floor(p*n + m))
        g = n*p + m - j
        if (g==0)&(j%2==0):
            gamma = 0.
        else:
            gamma = 1.
        return (1. - gamma)*x[j-1] + gamma*x[j]

    def four(x,p,n):
        m = 0.
        j = int(np.floor(p*n + m))
        gamma = n*p + m - j
        return (1. - gamma)*x[j-1] + gamma*x[j]

    def five(x,p,n):
        m = 0.5
        j = int(np.floor(p*n + m))
        gamma = n*p + m - j
        return (1. - gamma)*x[j-1] + gamma*x[j]

    def six(x,p,n):
        m = p
        j = int(np.floor(p*n + m))
        gamma = n*p + m - j
        return (1. - gamma)*x[j-1] + gamma*x[j]

    def seven(x,p,n):
        m = 1. - p
        j = int(np.floor(p*n + m))
        gamma = n*p + m - j
        return (1. - gamma)*x[j-1] + gamma*x[j]

    def eight(x,p,n):
        m = (p+1)/3.
        j = int(np.floor(p*n + m))
        gamma = n*p + m - j
        return (1. - gamma)*x[j-1] + gamma*x[j]

    def nine(x,p,n):
        m = p/4. + 3./8.
        j = int(np.floor(p*n + m))
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
    # Convert the percentile to a fraction
    if not fraction:
        p = p/100.
    if (p>1) or (p<0): raise InvalidPercentileError(p)
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

    def intrpl(y,j):
        if j>=n:
            return y[-1]
        elif j<1:
            return y[0]
        else:
            yj = y[j-1]
            yjp = y[j]
            return (1-f)*yj + f*yjp

    if y.ndim==1:
        y = np.sort(y)
        Qp = intrpl(y,j)
    else:
        spacedim = np.array(y.shape)[1:].prod()
        oldshape = y.shape[1:]
        y = y.reshape(y.shape[0], spacedim)
        Qp = np.ones(spacedim)*np.nan
        for lat in range(spacedim):
            array = np.sort(y[:,lat])
            Qp[lat] = intrpl(array,j)
        Qp = Qp.reshape(oldshape)
    return Qp


def quantile_zhang_fast(x, q, fraction=False, rmnans=False):
    """A faster version of quantile_zhang. Thanks goes to Mathias Hauser for this."""
    n = x.shape[0]
    # Convert the percentile to a fraction
    if fraction:
        q = q*100.
    if (q<0) or (q>100): raise InvalidPercentileError(q)
    if q==100.: return max(x)
    elif q==0.: return min(x)
    # adjust p such that np.percentile returns the same as quantile_zhang
    #  numpy uses (n -1); quantile_zhang uses (n + 1)
    q_adj = q * (n + 1.) / (n - 1.)
    #  numpy uses indices j and j+1; quantile_zhang uses j-1 and j
    q_adj -= 1. * 100./(n - 1)
    # ensure q_adj is in 0..100.
    q_adj = np.clip(q_adj, 0., 100.)

    return np.percentile(x, q_adj, axis=0)

def quantile_climpact(y,p,fraction=False):
    """quantile function used by climpact.

    Copy of the c_quantile function used in climpact.
    I have no idea where the interpolation is from.
    """
    if not fraction:
        p = p/100.
    if (p>1) or (p<0): raise InvalidPercentileError(p)
    y = np.array(y)
    if y.ndim==1:
        return qclimpact(y,p)
    else:
        spacedim = y.shape[1]
        nodims = y.ndim
        if nodims>2:
            spacedim = np.array(y.shape)[1:].prod()
            oldshape = y.shape[1:]
            y = y.reshape(y.shape[0], spacedim)
        Qp = np.ones(spacedim)*np.nan
        for lat in range(spacedim):
            if np.isnan(y[:,lat]).all(): continue
            Qp[lat] = qclimpact(y[:,lat],p)
        if nodims>2: Qp = Qp.reshape(oldshape)
    return Qp

def qclimpact(y,p):
    y = y[np.logical_not(np.isnan(y))]
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
    y = np.sort(y)
    if h==1:
        Qp = y[right_elem]
    elif h==0:
        Qp = y[left_elem]
    else:
        Qp = (1 - h)*y[left_elem] + h*y[right_elem]
    return Qp