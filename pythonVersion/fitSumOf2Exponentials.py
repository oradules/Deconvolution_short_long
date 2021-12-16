#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 12 11:43:43 2021

@author: rachel
"""
#    %%%%%%%%%%%%%%%%% 2-exp fit%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#    %%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
import numpy as np
from constraint0 import constraint0
from scipy.optimize import least_squares
def fitSumOf2Exponentials(xsgmin, xlmin, fsgmin, flmin, p1, p2, alpha):
    xsg=xsgmin
    xl=xlmin
    fsg=fsgmin
    fl=flmin
    alpha = alpha
    p2 = p2
    p1 = p1
    
    NS = len(xsg) 
    NL = len(xl)
    sNS=np.sqrt(NS)
    sNL=np.sqrt(NL) 
    fact1=np.sqrt(1-alpha)
    fact2=np.sqrt(alpha) 
    eshort=(1-fsg)*(1-p2)+p2;
    elong = (1-fl)*p1;
    
    if alpha == 1:
        def exp_fitness(k):
            return np.concatenate(([(k[2]*np.exp(k[0]*xsg)+(1-k[2])*np.exp(k[1]*xsg)-eshort)/sNS,(k[2]*np.exp(k[0]*xl)+(1-k[2])*np.exp(k[1]*xl)-elong)/sNL])) #%%%%% linear
    
    else:
        def exp_fitness(k):
            return np.concatenate((np.log((k[2]*np.exp(k[0]*xsg)+(1-k[2])*np.exp(k[1]*xsg))/eshort)/sNS*fact1,
                                           np.log((k[2]*np.exp(k[0]*xl)+ (1-k[2])*np.exp(k[1]*xl))/elong)/sNL*fact1,
                                           (k[2]*np.exp(k[0]*xsg)+(1-k[2])*np.exp(k[1]*xsg)-eshort)/sNS*fact2,
                                           (k[2]*np.exp(k[0]*xl)+ (1-k[2])*np.exp(k[1]*xl)-elong)/sNL*fact2)) ### mixed

    
    #%%%%%%% initial guess     
    k00= np.array([-0.1,-0.01,0.25], dtype= np.float64)
    k01= np.array([-0.1,-0.001,0.25])
    k0=np.zeros(3)
    amp = np.array([np.log(100),np.log(100)],dtype= np.float64) 
    
    NbIterationinFit=100    
    
    O2min=1e6    
    dofit=1  
    if dofit:
        for mc in range(NbIterationinFit):
            if mc%20 ==0:
                print(mc)
            k0=np.zeros(3)
            #### first try
            #### Change k00
            factor= np.exp(amp*(2*np.random.uniform(size=2)-1))
            k0[0:2] = k00[0:2]*factor
        
            #### sort ####
        
            k0[0:2]=np.sort(k0[0:2])
            k0min = k0[1]/(k0[1]-k0[0])
            k0max = 1
            k0[2] = k0min + np.random.uniform()*(k0max-k0min)
        
            e0 = exp_fitness(k0)
        
            if all(np.isfinite(e0)) and all(e0.imag==0):
            
                k = least_squares(exp_fitness, k0, ftol = (1e-8), max_nfev= 1e6, xtol= (1e-10)).x
                obj = sum(exp_fitness(k)**2)
            
                if obj < O2min and sum(k.imag)==0:
                    O2min=obj #### optimal objective function
                    kmin=k   #### optimal parameters

    return kmin, O2min