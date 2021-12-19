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
from fitSumOf2ExponentialsMain import fitSumOf2ExponentialsMain
from fitSumOf3ExponentialsContrained import fitSumOf3ExponentialsContrained
from scipy import interpolate
from scipy.optimize import least_squares

def fit2(dirwrite,name,dt,dtg,censored,censored_short,wt,wtc,lDEL,Total,Ninactive,visualize,time,sd):
    fsz=16 #figure size
    outliers_short = 1
    outliers_long = 1
    Ominmin=1e6
    optimalParams = []
    for ibig in range(3):
        alpha = round(0.9 - ibig * 0.3, 4)
        for cens in range(2):
            Kstore = np.empty((0,4)) # will store parameter and objective function values
            ######################## short movie ##################################
            if cens:
                xsg,fsg,logg,upg = ecdf_bounds(dtg,censored_short) #fsg will be the cumulative distribution function calculated at the points xsg using the data in dtg
            else:
                xsg,fsg,logg,upg = ecdf_bounds(dtg)
            #name = name + str(100*alpha)                   

            for nnn in range(0,7, 2):
                neg=0
                #print('nnn = {}'.format(nnn))
                wwt= np.append(wt,wtc)+3*(nnn)*60/2 ### shifted distribution from 0 t0 6 min        
                shift=3*nnn/2 #### shift in minutes        
                ########### estimate of integral for computing p1  #################
                a=3*60 # a is Pmin
                Tinactive =sum(wwt)
                Tactive = Total-Tinactive
                print('Tactive = ',Tactive, 'nnn = ', nnn)
                if Tactive<0:
                    nnn = nnn - 2
                    wwt= np.append(wt,wtc)+3*(nnn)*60/2 ### shifted distribution from 0 t0 6 min
                    shift=3*nnn/2 #### shift in minutes           
                    ########### estimate of integral for computing p1  #################
                    a=3*60 # a is Pmin
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
                print('estp1 ',estp1)
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
                #print(y1,y2)
                ##########################################################
                ################### Estimating p2  #######################
                p1=estp1 # this is p_l
                def objective(k):

                    return np.log((1-y2)*p1) - np.log((1-y1)*(1-k) + k)
                k0= 0.01
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
                    figfile=dirwrite+'/figure2_2exp_shift'+str(nnn)+'_cens'+str(cens)+'.pdf'
                    h.savefig(figfile)
                    plt.close()

                #%%%% fit distribution of spacings using combination of 2 exponentials. 3 
                #%%%% params


                kmin, OBJ, Kstore, optimalParams = fitSumOf2ExponentialsMain(xsg, xl, fsg, fl, p1, p2, alpha,cens, nnn, Kstore, Ominmin, optimalParams);
                #%optimalParams = {Ominmin ,kminmin,shiftmin,p1min,p2min,fsgmin,losmin,upsmin,flmin,lolmin,uplmin,xlmin,xsgmin,mw1opt,censmin}
                if optimalParams:
                    Ominmin=optimalParams[0];# %%%% optimal optimum
                    kminmin=optimalParams[1];# %%%% optimal optimal parameters
                    shiftmin=optimalParams[2];
                    p1min=optimalParams[3]; p2min=optimalParams[4];
                    fsgmin=optimalParams[5]; 
                    flmin=optimalParams[6];
                    xlmin=optimalParams[7]; 
                    xsgmin=optimalParams[8];
                    mw1opt=optimalParams[9];
                    alphamin=optimalParams[10];
                    censmin=optimalParams[11];
                    Kstoremin=optimalParams[12];
          
                l1=kmin[0];
                l2=kmin[1];
                A1=kmin[2];
                A2=1-A1;
                Omin = OBJ   
                if visualize:
                    h=plt.figure(105+nnn)
                    plt.loglog(xsg, (1-fsg)*(1-p2)+p2, label='short movie',marker='o',color='green',linestyle='', mfc='none')
                    plt.loglog(xl, (1-fl)*p1,label='long movie', marker='x', color='black',linestyle='')

                    sfit = kmin[2]*np.exp(kmin[0]*time)+(1-kmin[2])*np.exp(kmin[1]*time)
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
                    figfile=dirwrite+'/figure3_2exp_shift'+str(nnn)+'_cens'+str(cens)+'.pdf'
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
                sfit =  kmin[2]*np.exp(kmin[0]*time)+(1-kmin[2])*np.exp(kmin[1]*time)
                plt.loglog(time,sfit, label='fit 2 exp', color='black') ###fitted               
                plt.xlim([1e-1,1e5])
                plt.ylim([1e-6,1])
                plt.xlabel('\u0394 t [s]', fontsize=fsz)
                plt.ylabel('Survival function',fontsize=fsz)
                plt.legend()
                plt.title(label='cens='+ str(cens)+'\u0394\u2080='+str(int(shiftmin))+'Obj='+ str(round(Ominmin,5)),fontsize=fsz)
                plt.show
                figfile=dirwrite+'/figure3_optimal_2exp_'+name+'alpha'+str(alphamin)+'_cens'+str(censmin)+'.pdf'
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
                figfile=dirwrite+'/figure4_optimal_2exp_'+name+'_cens'+str(cens)+'.pdf'
                h.savefig(figfile)
                plt.close()
                fname = dirwrite+'/figure3_optimal_2exp_'+name+'alpha'+str(alpha)+'_cens'+str(int(cens))
                np.savez(fname+'.npz', xsgmin=xsgmin, fsgmin=fsgmin, xlmin=xlmin, flmin=flmin, time=time, sfit=sfit)    

            ###### optimal parameters    
            #%%%%%%%%%%%%%%%%%%%%%%%% optimal parameters %%%%%%%%%%%%%%%%%%%%%%%
            #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    kmin=kminmin;
    Kstore=Kstoremin
    l1=kmin[0];
    l2=kmin[1];
    A1=kmin[2];
    A2=1-A1;
    mean_w2 = -A1/l1-A2/l2; #%%%% average waiting time from parameters
    mean_w1=mw1opt;    

    S1=A1*l1+A2*l2;
    S2=A1*l1**2+A2*l2**2;
    S3=A1*l1**3+A2*l2**3;

    k2min=-S1;
    k1mmin=S1-S2/S1;
    k1pmin=(S3*S1-S2**2)/S1/(S1**2-S2);

    #%%%%%%% sort %%%%%%%%%%%%%%%%%%%%%%
    KK  = np.array([A1, A2]);
    K = np.array([l1,l2]);
    K_sorted = np.sort(K);
    isort = np.argsort(K)
    KK=KK[isort];
    l1=K_sorted[0];
    l2=K_sorted[1];
    A1=KK[0];
    A2=KK[1];
    l1min=l1;l2min=l2;A1min=A1;A2min=A2;

    res = np.array([k1pmin,k1mmin,k2min,l1min,l2min,A1min,A2min,Ominmin,45/mean_w2*60,shiftmin,alphamin,censmin,sd[1],sd[0]])

    #%%%%%%%%%%%%%%%% compute parameters with error bars %%%%%%%%%%%%%%
    ###################################################################
    overflow=1;# %%%%% 100% overflow    
    #%%%% select near-optimal parameters
    #%%%% O between Ominmin and Ominmin*(1+overflow)  %%%    
    Ksel = Kstore[ Kstore[:,3] < Ominmin*(1+overflow) , :  ]

    #%%%% compute intervals for parameters
    l1=Ksel[:,0];
    l2=Ksel[:,1];
    A1=Ksel[:,2];
    A2=1-A1;
    S1=A1*l1+A2*l2;
    S2=A1*l1**2+A2*l2**2;
    S3=A1*l1**3+A2*l2**3;
    K2=-S1;
    K1m=S1-S2/S1;
    K1p=(S3*S1-S2**2)/S1/(S1**2-S2);
    #          MRNA = 45/(-A1/l1-A2/l2-A3/l3)*60;


    #%%%%%%%%%%%%%%%%%%%%%%%% sort Ksel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Ksel_ext=np.transpose(np.array([A1,A2]));
    Ksel_sorted = np.sort(Ksel[:,0:2],1)
    isort = np.argsort(Ksel[:,0:2],1);
    nn=np.shape(isort);
    for i1 in range(nn[0]):
        Ksel_ext[i1,:] = Ksel_ext[i1,isort[i1,:]];
    Ksel_ext=np.hstack([Ksel_sorted,Ksel_ext]);
    l1=Ksel_ext[:,0];
    l2=Ksel_ext[:,1];
    A1=Ksel_ext[:,2];
    A2=Ksel_ext[:,3];
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    #######################################################################
    NS = len(xsgmin) 
    NL = len(xlmin)                

    resl= np.array([np.min(K1p),np.min(K1m),np.min(K2),np.min(l1),np.min(l2),np.min(A1),np.min(A2)])
    resl[0:5] = np.max(np.vstack([resl[0:5], np.zeros(5)]), axis=0)
    resh= np.array([np.max(K1p),np.max(K1m),np.max(K2),np.max(l1),np.max(l2),np.max(A1),np.max(A2)])

    return [res, resl, resh]
