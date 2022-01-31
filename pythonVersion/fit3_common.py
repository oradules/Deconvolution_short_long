import os
import numpy as np
import xlsxwriter
import shutil
from scipy.io import loadmat
from ecdf_estimate import ecdf_bounds
from scipy import interpolate
from scipy.optimize import least_squares
import matplotlib.pyplot as plt
import math
import random
from interpolateMetm import interpolateMetm
from fitSumOf3ExponentialsMainTest import *
#from fitSumOf3ExponentialsContrained import fitSumOf3ExponentialsContrained
from scipy import interpolate
from scipy.optimize import least_squares
from fit3_M1 import *
from fit3_M2 import *

def fit3_common(dirwrite,name,dt,dtg,censored,censored_short,wt,wtc,lDEL,Total,Ninactive,visualize,time,sd):
    print(sd)
    fsz=16 #figure size
    outliers_short = 1
    outliers_long = 1
    Ominmin=1e6
    optimalParams = []
    for ibig in np.arange(1):
        alpha = round(0.9 - ibig * 0.3, 4)
        for cens in range(2):
            Kstore = np.empty((0,6)) # will store parameter and objective function values
            ######################## short movie ##################################
            if cens:
                xsg,fsg,logg,upg = ecdf_bounds(dtg,censored_short) #fsg will be the cumulative distribution function calculated at the points xsg using the data in dtg
            else:
                xsg,fsg,logg,upg = ecdf_bounds(dtg)
            #name = name + str(100*alpha)                   

            for nnn in range(5):
                neg=0
                #print('nnn = {}'.format(nnn))
                wwt= np.append(wt,wtc)+3*(nnn)*60/2 ### shifted distribution from 0 t0 6 min        
                shift=3*nnn/2 #### shift in minutes        
                ########### estimate of integral for computing p1  #################
                a=1*60 # a is Pmin
                Tinactive =sum(wwt)
                Tactive = Total-Tinactive
                if Tactive<0:
                    nnn = nnn - 2
                    wwt= np.append(wt,wtc)+3*(nnn)*60/2 ### shifted distribution from 0 t0 6 min
                    shift=3*nnn/2 #### shift in minutes           
                    ########### estimate of integral for computing p1  #################
                    a=1*60 # a is Pmin
                    Tinactive =sum(wwt)
                    Tactive = Total-Tinactive
                    neg = 1

                Nactive = Ninactive - lDEL #### number of active periods
                imax=max(np.where(xsg<a)[0])
                imin= min(np.where(xsg>a)[0])

                x1=xsg[imin]
                f1=fsg[imin]
                x2=xsg[imax]
                f2=fsg[imax]
                f = f1 + (f2-f1)*(a-x1)/(x2-x1) # 1-S(a) based on the triangle function
                esp=(-a*(1-f)+np.trapz(1-fsg[xsg<a],xsg[xsg<a]))/f #espectation of the waiting time
                estp1= Ninactive/(Ninactive + Tactive/esp)

                ##################################################################
                ########### long movie

                if cens:
                    xl, fl, logg_l, upg_l=ecdf_bounds(wwt,censored)
                else:
                    xl, fl, logg_l, upg_l=ecdf_bounds(wwt)

                ############## fit p2 ###########################################
                fact1=1
                fact2=1
                uxsg, iu = np.unique(xsg, return_index=True)
                ufsg=fsg[iu]
                uxl, iu = np.unique(xl, return_index=True)
                ufl = fl[iu]
                M=max(uxsg)/fact2
                m=min(uxl[1:-1,])*fact1

                ##### interpolate on [m,M];
                y1,y2 = interpolateMetm(uxsg, ufsg, uxl, ufl, m, M, fact1, fact2)

                ##########################################################
                ################### Estimating p2  #######################
                p1=estp1 # this is p_l
                def objective(k):

                    return np.log((1-y2)*p1) - np.log((1-y1)*(1-k) + k)
                k0= 0.2
                k = least_squares(objective, k0, bounds=(0,1), ftol = (1e-15), max_nfev= 1e6, xtol= (1e-10)).x
                obj = sum(objective(k)**2)
                p2=k #this is ps
                #############################################################

                if visualize:                    
                    h=plt.figure(200+nnn)
                    plt.loglog(xsg, (1-fsg)*(1-p2)+p2, label='short movie', marker='o',color='green',linestyle='', mfc='none')
                    plt.loglog(xl, (1-fl)*p1,label='long movie',marker='o', color='black',linestyle='', mfc='none')
                    plt.xlabel('Time [s]', fontsize=fsz)
                    plt.ylabel('Frequency',fontsize=fsz)
                    plt.xlim([1e-1,1e4])
                    plt.ylim([1e-6,1e0])
                    plt.legend()
                    plt.title(label='\u0394\u2080='+str(3*nnn),fontsize=fsz)
                    plt.show
                    figfile=dirwrite+'/figure2_3exp_shift'+str(nnn)+'_cens'+str(cens)+'.pdf'
                    h.savefig(figfile)
                    plt.close()                  
                                        
                #%%%% fit distribution of spacings using combination of 2 exponentials. 3 
                #%%%% params

                kmin, OBJ, Kstore, optimalParams = fitSumOf3ExponentialsMainTest(xsg, xl, fsg, fl, logg, upg, logg_l,upg_l, p1, p2, alpha, nnn, Kstore, Ominmin, cens, optimalParams);
                #%optimalParams = {Ominmin ,kminmin,shiftmin,p1min,p2min,fsgmin,losmin,upsmin,flmin,lolmin,uplmin,xlmin,xsgmin,mw1opt,censmin}
                if optimalParams:
                    Ominmin=optimalParams[0];# %%%% optimal optimum
                    kminmin=optimalParams[1];# %%%% optimal optimal parameters
                    shiftmin=optimalParams[2];
                    p1min=optimalParams[3]; p2min=optimalParams[4];
                    fsgmin=optimalParams[5]; losmin=optimalParams[6];upsmin=optimalParams[7];
                    flmin=optimalParams[8]; lolmin=optimalParams[9];uplmin=optimalParams[10];
                    xlmin=optimalParams[11]; xsgmin=optimalParams[12];
                    mw1opt=optimalParams[13]; censmin=optimalParams[14];
                    alphamin=optimalParams[15];Kstoremin=optimalParams[16]
          
                l1=kmin[0];
                l2=kmin[1];
                A1=kmin[2];
                A2=kmin[3];
                A3=1-A1-A2;
                Omin = OBJ   
                if visualize:
                    h=plt.figure(105+nnn)
                    plt.loglog(xsg, (1-fsg)*(1-p2)+p2, label='short movie',marker='o',color='green',linestyle='', mfc='none')
                    plt.loglog(xl, (1-fl)*p1,label='long movie', marker='x', color='black',linestyle='')
                    sfit = kmin[3]*np.exp(kmin[0]*time)+kmin[4]*np.exp(kmin[1]*time)+(1-kmin[3]-kmin[4])*np.exp(kmin[2]*time)
                    
                    plt.loglog(time,sfit, label='fitted', color='red')    
                    #'mean_w2='+str("%.2g" % mean_w2)+' mRNA='+str("%.2g" % 45/mean_w2*60)
                    #plt.text(0.2,2e-4,'mean_w1='+str("%.2g" % mean_w1)
                    #plt.text(0.2,2e-5)                
                    plt.xlim([1e-1,1e5])
                    plt.ylim([1e-6,1])                
                    plt.xlabel('Time [s]', fontsize=fsz)
                    plt.ylabel('Frequency',fontsize=fsz)
                    plt.legend()
                    plt.title(label='cens='+ str(cens)+'\u0394\u2080='+str(int(shift))+'Obj='+ str(round(Omin,5)),fontsize=fsz)
                    plt.show
                    figfile=dirwrite+'/figure3_3exp_shift'+str(nnn)+'_cens'+str(cens)+'.pdf'
                    h.savefig(figfile)
                    plt.close()
                if neg==1:
                    break
            if visualize:
                ############ visualize optimal optimal fit 
                h=plt.figure(1000)
                plt.loglog(xsgmin, (1-fsgmin)*(1-p2min)+p2min, label='short movie',marker='o',color='red',linestyle='', mfc='none')
                plt.loglog(xlmin,(1-flmin)*p1min, label='long movie', marker='o', color='green',linestyle='',mfc='none')
                kmin=kminmin
                sfit = kmin[3]*np.exp(kmin[0]*time)+kmin[4]*np.exp(kmin[1]*time)+(1-kmin[3]-kmin[4])*np.exp(kmin[2]*time)
                plt.loglog(time,sfit, label='fit 3 exp', color='black') ###fitted               
                plt.xlim([1e-1,1e5])
                plt.ylim([1e-6,1])
                plt.xlabel('\u0394 t [s]', fontsize=fsz)
                plt.ylabel('Survival function',fontsize=fsz)
                plt.legend()
                plt.title(label='cens='+ str(cens)+'\u0394\u2080='+str(int(shiftmin))+'Obj='+ str(round(Ominmin,5)),fontsize=fsz)
                plt.show
                figfile=dirwrite+'/figure3_optimal_3exp_'+name+'alpha'+str(alphamin)+'_cens'+str(censmin)+'.pdf'
                h.savefig(figfile)
                plt.close()



                h=plt.figure(1003)
                plt.plot(xsgmin/60, (1-fsgmin)*(1-p2min)+p2min, marker='o',color='red',linestyle='', mfc='none')
                plt.plot(xlmin/60, (1-flmin)*p1min, marker='o',color='green',linestyle='', mfc='none')
                plt.loglog(time/60, sfit, color='black') ###fitted
            #            plt.xlim([0,3])
            #            plt.ylim([0,1])               
                plt.xlabel('\u0394 t [s]', fontsize=fsz)
                plt.ylabel('Survival function',fontsize=fsz)
                plt.legend()
                plt.title(label='cens='+ str(censmin)+'\u0394\u2080='+str(int(shiftmin))+'Obj='+ str(round(Ominmin,5)),fontsize=fsz)
                plt.show
                figfile=dirwrite+'/figure4_optimal_3exp_'+name+'alpha'+str(alphamin)+'_cens'+str(censmin)+'.pdf'
                h.savefig(figfile)
                plt.close()
                fname = dirwrite+'/figure3_optimal_3exp_'+name+'alpha'+str(alphamin)+'_cens'+str(int(censmin))
                np.savez(fname+'.npz', xsgmin=xsgmin, fsgmin=fsgmin, xlmin=xlmin, flmin=flmin, time=time, sfit=sfit)    

            ###### optimal parameters    
            #%%%%%%%%%%%%%%%%%%%%%%%% optimal parameters %%%%%%%%%%%%%%%%%%%%%%%
            #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [resM1, reslM1, reshM1] = fit3_M1(kminmin,Kstoremin,shiftmin,alphamin,censmin,sd)
    [resM2, reslM2, reshM2] = fit3_M2(kminmin,Kstoremin,shiftmin,alphamin,censmin,sd)
    return [resM1, reslM1, reshM1, resM2, reslM2, reshM2]