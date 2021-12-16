#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 12 18:01:06 2021

@author: rachel
"""

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%% constrained 3-exp fit %%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

import numpy as np
from scipy.optimize import least_squares
from constraint0 import constraint0
def fitSumOf3ExponentialsContrained(xsgmin, xlmin, fsgmin, flmin, p1, p2, alpha):
#### fit distribution of spacings using combination of 3 exponentials. 5 
    
    xsg=xsgmin;
    xl=xlmin;
    fsg=fsgmin;
    fl=flmin;
    NS = len(xsg); NL = len(xl);
    sNS=np.sqrt(NS);sNL=np.sqrt(NL); 
    eshort=(1-fsg)*(1-p2)+p2;
    elong = (1-fl)*p1;
    fact1 = np.sqrt(1-alpha);
    fact2 = np.sqrt(alpha);

    
    if alpha==1:
        def exp_fitness(k):
            if np.any(eshort==0):
                print('zero eshort')
                return
            if np.any(elong==0):
                print('zero elong')
                return
            return np.concatenate(([(k[3]*np.exp(k[0]*xsg)-(k[2]+k[3]*(k[0]-k[2]))/(k[1]-k[2])*np.exp(k[1]*xsg)+(k[1]+k[3]*(k[0]-k[1]))/(k[1]-k[2])*np.exp(k[2]*xsg)-eshort)/sNS,
                        (k[3]*np.exp(k[0]*xl)-(k[2]+k[3]*(k[0]-k[2]))/(k[1]-k[2])*np.exp(k[1]*xl)+(k[1]+k[3]*(k[0]-k[1]))/(k[1]-k[2])*np.exp(k[2]*xl)-elong)/sNL])) #%%%%% mixed
    else:
        def exp_fitness(k):
            if np.any(eshort==0):
                print('zero eshort')
            if np.any(elong==0):
                print('zero elong')
            return np.concatenate(([np.log((k[3]*np.exp(k[0]*xsg)-(k[2]+k[3]*(k[0]-k[2]))/(k[1]-k[2])*np.exp(k[1]*xsg)+(k[1]+k[3]*(k[0]-k[1]))/(k[1]-k[2])*np.exp(k[2]*xsg))/eshort)/sNS*fact1,
                                    np.log((k[3]*np.exp(k[0]*xl)-(k[2]+k[3]*(k[0]-k[2]))/(k[1]-k[2])*np.exp(k[1]*xl)+(k[1]+k[3]*(k[0]-k[1]))/(k[1]-k[2])*np.exp(k[2]*xl))/elong)/sNL*fact1,
                                    (k[3]*np.exp(k[0]*xsg)-(k[2]+k[3]*(k[0]-k[2]))/(k[1]-k[2])*np.exp(k[1]*xsg)+(k[1]+k[3]*(k[0]-k[1]))/(k[1]-k[2])*np.exp(k[2]*xsg)-eshort)/sNS*fact2,
                                    (k[3]*np.exp(k[0]*xl)-(k[2]+k[3]*(k[0]-k[2]))/(k[1]-k[2])*np.exp(k[1]*xl)+(k[1]+k[3]*(k[0]-k[1]))/(k[1]-k[2])*np.exp(k[2]*xl)-elong)/sNL*fact2])) #%%%%% mixed
        
    #%%%%%%% initial guess     
    k00= np.array([-0.1,-0.01,-0.001,0.25], dtype= np.float128)
    k0=np.zeros((4))
    amp = np.array([np.log(100), np.log(100), np.log(100)], dtype= np.float128)   
    NbIterationinFit=100      
    O3min=1e6    
    dofit=1  
    if dofit:
        test=0
        for mc in range(NbIterationinFit):
            if mc%20 ==0:
                print(mc)
            test=0
            while test==0: #test is just to re-do the iteration until we encounter no error
                #### first try
                #### Change k00
                factor= np.exp(amp*(2*np.random.uniform(size=3)-1))
                k0[0:3] = k00[0:3]*factor     
                #### sort ####
            
                if not ((k0[0] < k0[1]) and (k0[1] < k0[2])):
                    while not ((k0[0] < k0[1]) and (k0[1] < k0[2])):
                        factor = np.exp(np.multiply(amp,(2*np.random.rand(3)-1)))
                        k0[0:3] = np.multiply(k00[0:3],factor)
                    # end of while    
                #end of if
                k0[3] = 2*np.random.rand()-1
                if not( k0[1] + k0[3]*(k0[0]-k0[1]) < 0 ):
                    while not( k0[1] + k0[3]*(k0[0]-k0[1]) < 0 ):
                        k0[3]=2*np.random.rand()-1;  #%%%% try until condition satisfied
                
                e0 = exp_fitness(k0)
            
                if all(np.isfinite(e0)) and all(e0.imag==0):
                    try:
                        k = least_squares(exp_fitness, k0,bounds=(-np.inf,[0,0,0,np.inf]),ftol = (1e-8),max_nfev= 1e6, xtol= (1e-10)).x
                        obj = sum(exp_fitness(k)**2)
                        test=1
                    except:
                        pass        
                    if obj < O3min and sum(k.imag)==0:
                        O3min=obj #### optimal objective function
                        kmin=k   #### optimal parameters

    return kmin, O3min