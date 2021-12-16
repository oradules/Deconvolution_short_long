#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 11 21:36:32 2021

@author: rachel
"""
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
def  interpolateMetm(uxsg, ufsg, uxl, ufl, m, M, fact1, fact2):
    #%%%%% interpolate on [m,M];
    xss = uxsg
    yss = ufsg
    src = [0,0]; #% (0,0)
    dest = [max(m,M),1]; #% (m,1)
    
    xss = np.hstack((xss,dest[0])) 
    yss = np.hstack((yss,dest[1]))
    t = np.linspace(0, 1, np.size(xss))
    
    mn=5000
    tt = np.linspace(0, 1, mn);
    xx = interpolate.interp1d(t,xss, kind = 'cubic')(tt)
    yy = interpolate.interp1d(t,yss, kind = 'cubic')(tt)
    

    # plt.plot(xss,yss,'ro');
    # plt.plot(xx, yy, 'b', 'LineWidth', 1.5);
    
    xll = uxl
    yll = ufl;
    t = np.linspace(0, 1, np.size(xll))
    
    mn=5000
    ttl = np.linspace(0, 1, mn);
    xxl =  interpolate.interp1d(t,xll, kind = 'cubic')(ttl)
    yyl =  interpolate.interp1d(t,yll, kind = 'cubic')(ttl)
    

    # plt.plot(xss,yss,'ro')
    # plt.plot(xxl, yyl, 'b', 'LineWidth', 1.5);
    # plt.plot(xll,yll, 'ro');
    # plt.plot(xx, yy, 'g', 'LineWidth', 1.5);
    
    
    M_ = max(uxsg)/fact2;
    m_ = min(uxl[1:])*fact1;
    m = min(m_,M_); M = max(m_,M_);
    inddm = np.transpose(np.logical_and([xxl>m],[xxl<M])[0])
    
    y2 = yyl[inddm];xt = xxl[inddm]
    
    indxx1 = np.where(xx>=m)[0];
    indxx1 = indxx1[0]-1;
    indxx2 = np.where(xx<=M)[0];
    indxx2 = indxx2[-1];
    ytest = yy[indxx1:indxx2+1];
    x1 = xx[indxx1:indxx2+1];
    
    #[ux1,ux1i] = unique(x1);
    ux1, ux1i = np.unique(x1, return_index=True)
    x1 = x1[ux1i];
    ytest = ytest[ux1i];
    y1 = interpolate.interp1d(x1, ytest, kind = 'cubic')(xt)

    return y1, y2