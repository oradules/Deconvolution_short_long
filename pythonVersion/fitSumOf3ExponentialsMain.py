#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 12 14:04:23 2021

@author: rachel
"""
#    %%%%%%%%%%%%%%%%% 3-exp fit%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#    %%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
import numpy as np
from scipy.optimize import least_squares
from constraint0 import constraint0
def fitSumOf3ExponentialsMain(xsg, xl, fsg, fl, los, ups, lol, upl, p1, p2, alpha, nnn, Kstore, Ominmin, cens, optimalParams):
#### fit distribution of spacings using combination of 3 exponentials. 5 
    
    xl=xl[1:-1]
    fl=fl[1:-1]

    xsg = xsg[1:-1]
    fsg = fsg[1:-1]
    
    los = los[1:-1]
    ups = ups[1:-1]
    lol=lol[1:-1];
    upl=upl[1:-1];
    
    NS = len(xsg) 
    NL = len(xl)
    sNS=np.sqrt(NS)
    sNL=np.sqrt(NL) 
    fact1=np.sqrt(1-alpha)
    fact2=np.sqrt(alpha) 
    eshort=(1-fsg)*(1-p2)+p2;
    elong = (1-fl)*p1;

#%%%%%%% compute average waiting time from AUC of the distribution
    mean_w1=np.trapz(np.append(eshort,elong),np.append(xsg,xl))
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if alpha == 1: #%%% linear 
        def exp_fitness(k):
            return np.concatenate(([(k[3]*np.exp(k[0]*xsg)+k[4]*np.exp(k[1]*xsg)+(1-k[3]-k[4])*np.exp(k[2]*xsg)-eshort)/sNS,
                                    (k[3]*np.exp(k[0]*xl)+k[4]*np.exp(k[1]*xl)+(1-k[3]-k[4])*np.exp(k[2]*xl)-elong)/sNL])) #; %%%%% linear
    else:
        def exp_fitness(k):
            ee = np.concatenate(((np.log(((k[3]*np.exp(k[0]*xsg)+k[4]*np.exp(k[1]*xsg)+(1-k[3]-k[4])*np.exp(k[2]*xsg))/eshort))/sNS*fact1,
                                               np.log((k[3]*np.exp(k[0]*xl)+k[4]*np.exp(k[1]*xl)+(1-k[3]-k[4])*np.exp(k[2]*xl))/elong)/sNL*fact1,
                                               ((k[3]*np.exp(k[0]*xsg)+k[4]*np.exp(k[1]*xsg)+(1-k[3]-k[4])*np.exp(k[2]*xsg))-eshort)/sNS*fact2,
                                               ((k[3]*np.exp(k[0]*xl)+k[4]*np.exp(k[1]*xl)+(1-k[3]-k[4])*np.exp(k[2]*xl))-elong)/sNL*fact2))) ### mixed
            return ee
 
    
    #%%%%%%% initial guess     
    k00= np.array([-0.1,-0.01,-0.001,0.25,0.25], dtype= np.float64)
    k0=np.zeros((5))
    amp = np.array([np.log(100),np.log(100), np.log(100)], dtype= np.float32)   
    NbIterationinFit=100      
    Omin=1e6    
    dofit=1  
    if dofit:
        for mc in range(NbIterationinFit):
            test=0
            print(mc)
            while test==0: #test is just to re-do the iteration until we encounter no error
                print('trying least squares for mc = {}'.format(mc))
                #### first try #### Change k00
                factor= np.exp(amp*(2*np.random.uniform(size=3)-1))
                k0[0:3] = k00[0:3]*factor           
                #### sort ####
                if not ((k0[0] < k0[1]) and (k0[1] < k0[2])):
                    while not ((k0[0] < k0[1]) and (k0[1] < k0[2]) and (not np.any(np.isinf(k0[0:3])))):
                        factor = np.exp(amp*(2*np.random.rand(3)-1))
                        k0[0:3] = k00[0:3]*factor
                        print('while 1')
                        #print(k0)
                    # end of while    
                #end of if
                k0[3:] = 2*np.random.rand(2)-1    
                if not constraint0(k0):
                    while not constraint0(k0):
                        k0[3:] =  2*np.random.rand(2)-1 
                        print('while 2')
                e0 = exp_fitness(k0)           
                if all(np.isfinite(e0)) and all(e0.imag==0):
                    print('sorting')
                    Atemp = np.asarray([k0[3],k0[4], 1-k0[3]-k0[4]]) # A values before sorting 
                    kk_temp = np.sort(k0[0:3])
                    IX = np.argsort(k0[0:3])
                    k0[0:3]= kk_temp
                    Atemp=Atemp[IX]
                    k0[3:5]=Atemp[0:2]                    
                    try:
                        print('trying lsq')
                        k = least_squares(exp_fitness, k0, bounds=(-np.inf,[0,0,0,np.inf,np.inf]), ftol = (1e-8), max_nfev= 1e6, xtol= (1e-10)).x
                        obj = sum(exp_fitness(k)**2)
                        test = 1
                        print('value found!')                        
                    except:
                        pass
                    Atemp = np.asarray([k[3],k[4], 1-k[3]-k[4]]) # A values before sorting 
                    kk_temp = np.sort(k[0:3])
                    IX = np.argsort(k[0:3])
                    k[0:3]= kk_temp
                    Atemp=Atemp[IX]
                    k[3:5]=Atemp[0:2]
            
                    if obj < Omin and sum(k.imag)==0:
                        Omin=obj #### optimal objective function
                        kmin=k   #### optimal parameters
                        shift=3*nnn/2 #### shift in minutes
                        p1min=p1 #### p1 value
                        p2min=p2 #### p2 value
                    if (sum(k.imag)==0) and constraint0(k):
                        Kstore=np.vstack((Kstore,np.concatenate((k,[obj])))) #%%%% consider suboptimal values

        if Omin < Ominmin:
            Ominmin=Omin;# %%%% optimal optimum
            kminmin=kmin;# %%%% optimal optimal parameters
            shiftmin=shift;
            p1min=p1;
            p2min=p2;
            fsgmin=fsg; losmin=los; upsmin=ups;
            flmin=fl; lolmin=lol; uplmin=upl;
            xlmin=xl;
            xsgmin=xsg;
            mw1opt=mean_w1;
            censmin=cens;
            alphamin = alpha;
            Kstoremin = Kstore;
            optimalParams = [Ominmin ,kminmin,shiftmin,p1min,p2min,fsgmin,losmin,upsmin,flmin,lolmin,uplmin,xlmin,xsgmin,mw1opt,censmin,alphamin,Kstoremin];   
        OBJ = Omin 
    return  kmin, OBJ, Kstore, optimalParams
