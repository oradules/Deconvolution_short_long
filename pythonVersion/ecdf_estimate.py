#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  8 14:57:54 2020

Function to estimate the pointwise 
95% lower and upper confidence bounds evaluated at the points x 
using the censored data 'y' specified by 'censored'

y : array of values
censored : array of 0 and 1 where 1 indicates censored value.
Returns: Lower_bounds, Upper_bounds

By default alpha is set to 0.05 to calculate 95% confidence bound


@author: Rachel Topno
"""

import numpy as np
import math
from scipy.stats import norm

def ecdf_bounds(y, censored=np.asarray([])):
    fn = 'cdf'
    talpha = 0.05
    cdf_sf = 1
    ty = y
    tx = np.asarray(ty)
    censored
    if censored.size == 0:
        tcen = np.zeros(tx.size)
    else:
        tcen = censored


    tn = len(tx)
    tt = np.argsort(tx)
    tx = tx[tt, ]
    tcen = tcen[tt, ]
    tfreq = np.ones(tx.size)
    totcumfreq = np.cumsum(tfreq)  # how many waiting times there are
    # np.cumsum(np.multiply(tfreq, np.multiply(np.logical_not(tcen),1))) # how many uncencored waiting times
    ob_cumfreq = np.cumsum(np.multiply(tfreq, np.multiply(np.logical_not(tcen), 1)))  # how many uncencored witing times
    tt = np.where(np.diff(tx) == 0)[0]  # find indices where value is changing
    if np.any(tt):
        tx = np.delete(tx, tt)  # %% remove values which are repeating consecutively
        totcumfreq = np.delete(totcumfreq, tt)
        ob_cumfreq = np.delete(ob_cumfreq, tt)

    totcount = totcumfreq[-1]  # total waiting times


    tD = np.concatenate(([ob_cumfreq[0]], np.diff(ob_cumfreq)))  # time to death (activation)
    tN = totcount - np.concatenate(([0], totcumfreq[0:-1]))  # at risk times
    tt = np.where(tD > 0)[0]  # remove obs where time to death is zero
    tx = tx[tt, ]
    tD = tD[tt, ]
    tN = tN[tt, ]

    if cdf_sf:  # % 'cdf' or 'survivor'
        tS = np.cumprod(1 - np.divide(tD, tN))
        if fn in 'cdf':
            tFunc = 1 - tS
            tF0 = 0
        else:
            tFunc = np.cumsum(np.divide(tD, tN))
            tF0 = 1
    else:
        tFunc = np.cumsum(np.divide(tD, tN))
        tF0 = 0

    tx = np.concatenate(([min(ty)], tx))
    tF = np.concatenate(([tF0], tFunc))

    # % Get standard error of requested function
    if cdf_sf:  # 'cdf' or 'survivor'
        tse = np.empty(tD.size)
        tse[:] = np.NaN
        if tN[-1] == tD[-1]:
            tt = np.arange(0, len(tN) - 1)
        else:
            tt = np.arange(0, len(tN))
        tse[tt] = np.multiply(tS[tt], np.sqrt(np.cumsum(np.divide(tD[tt], (np.multiply(tN[tt], (tN[tt]-tD[tt])))))))
    else:  # % 'cumhazard'
        tse = np.sqrt(np.cumsum(np.divide(tD, np.multiply(tN, tN))))

    # % Get confidence limits
    zalpha = -norm.ppf(talpha / 2)
    halfwidth = zalpha*tse
    Flo = np.maximum([0], tFunc-halfwidth)

    # Flo(np.isnan(halfwidth)) = NaN; % max drops NaNs, put them back
    if cdf_sf:  # % 'cdf' or 'survivor'
        Fup = np.minimum(1, tFunc+halfwidth)

    #    Fup(isnan(halfwidth)) = NaN; % max drops NaNs
    else:  # % 'cumhazard'
        Fup = tFunc + halfwidth  # % no restriction on upper limit
    Flo = np.concatenate(([math.nan], Flo))
    Fup = np.concatenate(([math.nan], Fup))


    return tx, tF, Flo, Fup

