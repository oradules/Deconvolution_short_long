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
def fitSumOf2ExponentialsMain(xsg, xl, fsg, fl, p1, p2, alpha, cens, nnn, Kstore, Ominmin, optimalParams):
    print('Fitting sum of two exponentials...') 
    xl=xl[1:-1]
    fl=fl[1:-1]

    xsg = xsg[1:-1]
    fsg = fsg[1:-1]
    
    # los = los[1:-1]
    # ups = ups[1:-1]
    # lol=lol[1:-1];
    # upl=upl[1:-1];
    
    NS = len(xsg) 
    NL = len(xl)
    sNS=np.sqrt(NS)
    sNL=np.sqrt(NL) 
    fact1=np.sqrt(1-alpha)
    fact2=np.sqrt(alpha) 
    eshort=(1-fsg)*(1-p2)+p2;
    elong = (1-fl)*p1;
    
    mean_w1=np.trapz(np.append(eshort,elong),np.append(xsg,xl))    
    
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
    k0=np.zeros(3)
    amp = np.array([np.log(100),np.log(100)],dtype= np.float64) 
    
    NbIterationinFit=100    
    
    Omin=1e6    
    dofit=1  
    if dofit:
        for mc in range(NbIterationinFit):
            test=0
            while test==0: #test is just to re-do the iteration until we encounter no error
                #print('trying least squares for mc = {}'.format(mc))
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
                    try:
                        k = least_squares(exp_fitness, k0, ftol = (1e-8), max_nfev= 1e6, xtol= (1e-10)).x
                        obj = sum(exp_fitness(k)**2)
                        test = 1
                    except:
                        pass
                    if (obj < Omin) and sum(k.imag)==0:
                       Omin=obj; #%%%% optimal objective function
                       kmin=k;   #%%%% optimal parameters
                       shift=3*nnn/2; #%%%% shift in minutes
                       p1min=p1; #%%%% p1 value
                       p2min=p2; #%%%% p2 value
                    if sum(k.imag)==0:
                        Kstore = np.vstack((Kstore,np.concatenate((k,[obj]))))#%%%% consider suboptimal values

    if Omin < Ominmin:
           Ominmin=Omin; #%%%% optimal optimum
           kminmin=kmin; #%%%% optimal optimal parameters
           p1min=p1;
           p2min=p2;
           fsgmin=fsg;
           flmin=fl;
           xlmin=xl;
           xsgmin=xsg;
           mw1opt=mean_w1;
           shiftmin = shift;
           alphamin = alpha;
           censmin = cens;
           Kstoremin = Kstore;
           optimalParams = [Ominmin ,kminmin,shiftmin,p1min,p2min,fsgmin,flmin,xlmin,xsgmin,mw1opt,alphamin,censmin,Kstoremin];
    return kmin, Omin, Kstore, optimalParams
