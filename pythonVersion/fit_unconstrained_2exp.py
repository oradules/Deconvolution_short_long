#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 11 16:20:35 2021

@author: rachel
"""
import sys
sys.path.append('/home/rachel/Documents/longMovie/long_movie_artificial/')
#########################################################################
##### Code written by Ovidiu Radulescu, University of Montpellier, June 2019
##### reads results of genetic algorithm for short movies, data from long
##### movies
##### needs 1) short movies decomvolution results  result_tat_name_short.mat
##### 2) long movie data name_long.mat or long movie data name_long_raw.mat
##### performs unconstrained two exponentials fit, compute parameters of
##### the two states model ##############################################
#########################################################################

##### Library needed ########################################################
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
from constraint0 import constraint0
from interpolateMetm import interpolateMetm
from fitSumOf2ExponentialsMain import fitSumOf2ExponentialsMain
from fitSumOf3Exponentials import fitSumOf3Exponentials
from fitSumOf3ExponentialsContrained import fitSumOf3ExponentialsContrained
from posPred2cPosPred import *
########################################################################

os.chdir('/home/rachel/Documents/longMovie/long_movie_artificial/') 
os.listdir()
fsz=16 ### font size for the figures
outliers_long=1 ##### this is to change directly if we dont want to comput the outliners for long movies
msz=5
lw=1
outliers_short=0 ####  this is to change directly if we dont want to comput the outliners for short movies
visualize = 1 #### if we dont want to draw the graph
#fact=0.75 %%% [fact/(1-fact)]^2 is the relative importance of linear scale in the objective function

##### file names ########################################################
#%%%% file names %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
list_short = ['result_HIVlo.npz'];
files_long = ['HIVlo_long.mat'];
names = ['HIVlo'];
# list_short = ['result_clean_eef1a1_wtcPosPred10.mat','result_clean_eef1a1_tatamutcPosPred10.mat','result_clean_eef1a1_g4mutcPosPred10.mat']
# files_long  = ['data_long_movies_eef1a1_wt.mat','data_long_movies_eef1a1_tatamut.mat','data_long_movies_eef1a1_g4mut.mat']
# names = ['eef1a1_wt','eef1a1_tatamut','eef1a1_g4mut']



############# Parameters ###############################################
TaillePreMarq = 700 # 700 bases (ref publi)
TailleSeqMarq = 2900 # bases (ref publi)
EspaceInterPolyMin = 30 # space between polymerase in bases (orientative)
Polym_speed = 67 # average speed bases per second (Ref publi)
TaillePostMarq = 1600+100*Polym_speed # 1600 bases (ref publi)
FreqEchImg = (1/3) # 1/3 image per second data time sampling
Intensity_for_1_Polym = 1 ###
DureeSignal = (TaillePreMarq + TailleSeqMarq + TaillePostMarq) / Polym_speed # signal length in seconds

#### other parameters short movie ######################################

lframes=400 ## number of frames
tstep = 3 ### frames period in seconds
DureeSimu = lframes*tstep ### movie length in seconds
frame_num = lframes ### number of frames
DureeAnalysee = DureeSignal + DureeSimu #bcz the signal may continue after the simulation is done
FreqEchSimu = 1/(EspaceInterPolyMin/Polym_speed) # how many intervals (possible poly start position) which is equivalent to max number of polymerase 
num_possible_poly = round(DureeAnalysee/(EspaceInterPolyMin/Polym_speed)) #

########################################################################

tmax=DureeAnalysee ### duration of timeseries (consider all the polymerases between 0 and tmax)

########################################################################

for ibig in range(3):
    alpha = round(0.9 - ibig * 0.3, 4)


### this part is to create the folder where we want to put the result and delete the old folder before when we want to run it again 
    DataFilePath0 = 'resultFit2Exp'#+str(round(100*alpha)) #### where to write results
    if os.path.exists(DataFilePath0):
        
        for filename in DataFilePath0:
            file_path = os.path.join(DataFilePath0, filename)
            try:
                if os.path.isfile(file_path) or os.path.islink(file_path):
                    os.unlink(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)
            except Exception as e:
                print('Failed to delete %s. Reason: %s') # (file_path, e))
          
    else:
        os.mkdir(DataFilePath0)
    
############# creating an xls file to put the data in it ###############################    
    xlsfilename = DataFilePath0 + '/results_2exp_' + str(100*alpha)+ '.xlsx'    
    workbook = xlsxwriter.Workbook(xlsfilename) # Workbook() takes one, non-optional, argument which is the filename that we want to create.
    worksheet = workbook.add_worksheet() # The workbook object is then used to add new worksheet via the add_worksheet() method.
    
    # 'Fname','OBJ','shift','cens','lambda1','lambda2','A1','A2','k2','k1m','k1p','mRNA'
    
    # Use the worksheet object to write data via the write() method
    worksheet.write('A1', 'Fname') 
    worksheet.write('B1', 'OBJ')    
    worksheet.write('C1', 'shift') 
    worksheet.write('D1', 'cens')
    worksheet.write('E1', 'lambda1')
    worksheet.write('F1', 'lambda2')
    worksheet.write('G1', 'A1')
    worksheet.write('H1', 'A2')
    worksheet.write('I1', 'K1p')
    worksheet.write('J1', 'K1m')
    worksheet.write('K1', 'K2')
    worksheet.write('L1', 'mRNA')

    ################################################################
    for iname in range(len(names)):
        name=names[iname]
        files_short=list_short[iname]
        fnameL =files_long[iname] 
    
        ######## dirwrite is a folder inside the folder that we already created so that we dont put the excel sheet with the images as this will contain the images
        dirwrite= DataFilePath0+'/'+name+'_result'+str(round(100*alpha)) #### where to write result

        if os.path.exists(dirwrite):
            shutil.rmtree(dirwrite, ignore_errors = True)

        os.mkdir(dirwrite)

    ################ loading data obtained from the ga for the short movies 
    ################ and the data that we have from the long movie  and 
    ################ concatenate if there are severale file #####
    ########### for short movies ########################

        fnameS = files_short
        load_dataS =np.load(fnameS, allow_pickle=True)   #### load datafile
        
        DataExp = load_dataS['DataExp']
        DataPred= load_dataS['DataPred']
        PosPred= load_dataS['PosPred']
        Fit= load_dataS['Fit']
#        load_traj=loadmat('lowtatartificialbis_positions.mat') 
        cPosPred = np.asarray(posPred2cPosPred(PosPred))

        #fnameL = files_long
        
        load_dataL =loadmat(fnameL)   #### load datafile
        DataExpLong = load_dataL['DataExpLong']
        DataExpLongRaw = load_dataL['DataExpLongRaw']
        Time = load_dataL['Time']
        Time = Time.reshape(np.size(Time),)
        #Time = np.arange(0, len(DataExpLong)*3,3)
    
        n_cells = len(DataPred[0]) #### number of cells for each cell_type (mut or tatmut...)

######### eliminate outliers handling short movie ########


#""""Several observed transcription sites had abnormal behaviour (too many or too few events). The decision
# was made to take them off the data set as follows which is what we are doing with outliners handling"""" 

        if outliers_short: 

            Ev =np.zeros((n_cells)) ####  based on quartile of the numbers of polymerases which is equivalent to saying: Compute the amount of events of transcription that happened during the movie
            for i in range(n_cells): #### for each cell
                Ev[i] = len(cPosPred[0,i][0]) #### estimate number of polymerase for each cell
            Q = np.quantile(Ev, [0.25, 0.5, 0.75]) ### devide Ev into 3 quantile

            isel = np.where((Q[0]-1.5*(Q[2]-Q[0]) < Ev) & (Ev < Q[2] + 1.5*(Q[2]-Q[0])))[0]  ####we are restricting our computation to this interval
        
            DataExp=DataExp[:,isel]
            DataPred=DataPred[:,isel]
            cPosPred=cPosPred[:,isel]
            #Fit=Fit[isel]
            n_cells=len(isel)

##########################################################


######## compute distribution of spacings short movie ######

        dt = np.array([])
        dtc = np.array([]) #"""" dtc represent what """
        tmax=0

        for i in range (n_cells):
            indices = cPosPred[0,i][0] #### indice is the estimate position of polymerase for each cell
            times = indices / FreqEchSimu ### last transcription in the movie
            Mtimes=max(times)
        
            if Mtimes > tmax:
                tmax = Mtimes

        for i in range (n_cells):
            indices = cPosPred[0,i][0]
            times = indices / FreqEchSimu
            
            lt = len(times)
            
            if lt == 0:
                dtc = np.append(dtc,tmax)
            elif lt ==1:
                dtc = np.append(dtc,(tmax-times[0]))
            else:
                dtimes = np.diff(times) #### uncensured intervals
                dt = np.append(dt, dtimes)
            #### censored intervals
                if tmax>times[-1]: #### cz we dont know the last waiting time for times[-1]
                    dtc= np.append(dtc,(tmax-times[-1]))
                if times[0]>0: ### bcz we don't know how much times has passed this the polymerase that occured before the begining of the movie and the first one in the movie
                    dtc = np.append(dtc,times[0])
                
                
#### eliminate zero intervals #########
        dt = np.delete(dt, np.where(dt==0))

        dtg= np.append(dt,dtc)

        censored_short=np.append(np.zeros((len(dt))), np.ones((len(dtc))))


###### eliminate outliers handling long movie #######
        '''
"""""""""" Some transcription sites in long movies data set also show unusual behaviours being active (FI=0) or inactive
(FI=1) during the entire movie. We exclude these outliers as we did it for the short movies, but based on the
fraction of inactivity for each transcription site.""""""""""""""""
        '''
        if outliers_long:
            FI=np.zeros(len(DataExpLong[0])) ## fraction of inactivity
            for i in range(len(DataExpLong[0])):
            
                N0=len(np.where(DataExpLong[:,i]==0)[0]) ### number of indices where we don't have a polymerase
                N1=len(np.where(DataExpLong[:,i]==1)[0]) ### number of indices where we have polymerase
            
                FI[i] = N0/(N0+N1) 

            Q = np.quantile(FI, [0.25, 0.5, 0.75])
            isel = np.where((Q[0]-1.5*(Q[2]-Q[0]) < FI) & (FI < min(Q[2] + 1.5*(Q[2]-Q[0]),1)) & (FI>0))[0]
            DataExpLong=DataExpLong[:,isel] ####we are restricting our computation to this interval

########################################################################

        store = DataExpLong

###### estimate distributions for long movies ##########################
        wt=np.array([])
        wtc=np.array([])
        Ninactive=0
        Total=0
        tmaxlomax=0

        for i in range(len(store[0])): ### all cell

            ilast = np.where(store[:,i]==2)[0]
        
            if ilast.size == 0:
                tmaxlo=Time[-1]*60
            else:
                tmaxlo=Time(min(ilast))*60

            ind = np.where(store[:,i]==1)[0]

            if ind.size !=0:
                laps = np.diff(ind)
                find_laps=np.where(laps>1)[0]
                Ninactive = Ninactive+len(find_laps)+1
                wtimes=(laps[laps>1]-1)*tstep*60
                wlast=tmaxlo - (1+ind[-1])*tstep*60 #### tmaxlo depends on cell
                wt=np.append(wt, wtimes) ### waiting times in seconds

                if wlast>0:
                    wtc = np.append(wtc, wlast) ### add last interval
            Total = Total + tmaxlo
            if tmaxlo > tmaxlomax:
                tmaxlomax = tmaxlo

        time = np.arange(start=0, stop=tmaxlomax+1, step=0.1) ### time in seconds
        censored=np.append(np.zeros(len(wt)),np.ones(len(wtc))) ### censored long movie
        
        for cens in range(2):
            
            Ominmin=1e6
    
            Kstore = np.empty((0,4)) # will store parameter and objective function values
            optimalParams = []
    ######################## short movie ##################################
    
            if cens:
                xsg,fsg,logg,upg =ecdf_bounds(dtg,censored_short) #fsg will be the cumulative distribution function calculated at the points xsg using the data in dtg
            else:
                xsg,fsg,logg,upg = ecdf_bounds(dtg)
    
            name = name + str(100*alpha)
            
            for nnn in range(0,7, 2):
                print('nnn = {}'.format(nnn))
                wwt= np.append(wt,wtc)+3*(nnn)*60/2 ### shifted distribution from 0 t0 6 min        
                shift=3*nnn/2 #### shift in minutes        
                ########### estimate of integral for computing p1  #################
                a=3*60 # a is Pmin
                Tinactive =sum(wwt)
                Tactive = Total-Tinactive
                if Tactive<0:
                    nnn = nnn - 2
                    wwt= np.append(wt,wtc)+3*(nnn)*60/2 ### shifted distribution from 0 t0 6 min
                    shift=3*nnn/2 #### shift in minutes           
                    ########### estimate of integral for computing p1  #################
                    a=3*60 # a is Pmin
                    Tinactive =sum(wwt)
                    Tactive = Total-Tinactive
                    
                Nactive = Ninactive - len(DataExpLong[0]) #### number of active periods
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
                    figfile=dirwrite+'/figure2_2exp_shift'+str(nnn)+'_cens'+str(cens)+'.pdf'
                    h.savefig(figfile)
                    plt.close()
                    
                #%%%% fit distribution of spacings using combination of 3 exponentials. 5 
                #%%%% params
      
                kmin, OBJ, Kstore, optimalParams = fitSumOf2ExponentialsMain(xsg, xl, fsg, fl, p1, p2, alpha, nnn, Kstore, Ominmin, optimalParams);
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
    

            if visualize:
                ############ visualize optimal optimal fit 
                h=plt.figure(1000)
                plt.loglog(xsgmin, (1-fsgmin)*(1-p2min)+p2min, label='short movie',marker='o',color='red',linestyle='', mfc='none')
                plt.loglog(xlmin,(1-flmin)*p1min, label='long movie', marker='o', color='green',linestyle='',mfc='none')
                kmin=kminmin
                sfit =  kmin[2]*np.exp(kmin[0]*time)+(1-kmin[2])*np.exp(kmin[1]*time)
                plt.loglog(time,sfit, label='fit 3 exp', color='black') ###fitted               
                plt.xlim([1e-1,1e5])
                plt.ylim([1e-6,1])
                plt.xlabel('\u0394 t [s]', fontsize=fsz)
                plt.ylabel('Survival function',fontsize=fsz)
                plt.legend()
                plt.title(label='cens='+ str(cens)+'\u0394\u2080='+str(int(shift))+'Obj='+ str(round(Omin,5)),fontsize=fsz)
                plt.show
                figfile=dirwrite+'/figure3_optimal_2exp_'+name+'_cens'+str(cens)+'.pdf'
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
                plt.title(label='cens='+ str(cens)+'\u0394\u2080='+str(int(shift))+'Obj='+ str(round(Omin,5)),fontsize=fsz)
                plt.show
                figfile=dirwrite+'/figure4_optimal_2exp_'+name+'_cens'+str(cens)+'.pdf'
                h.savefig(figfile)
            
                plt.close
                fname = dirwrite+'/figure3_optimal_2exp_'+name+'_cens'+str(int(cens))
                np.savez(fname+'.npz', xsgmin=xsgmin, fsgmin=fsgmin, xlmin=xlmin, flmin=flmin, time=time, sfit=sfit)    
            
            ###### optimal parameters    
            #%%%%%%%%%%%%%%%%%%%%%%%% optimal parameters %%%%%%%%%%%%%%%%%%%%%%%
            #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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
    
            #%%%%%%%%%%%%%%%% compute parameters with error bars %%%%%%%%%%%%%%
            ###################################################################
            overflow=1;# %%%%% 100% overflow    
            #%%%% select near-optimal parameters
            #%%%% O between Ominmin and Ominmin*(1+overflow)  %%%    
            Ksel = Kstore[ Kstore[:,3] < Ominmin*(1+overflow) , :  ]
            #%%%% compute intervals for parameters

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
    
 #           MRNA = 45/(-A1/l1-A2/l2-A3/l3)*60;
    

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
            
            xlsfilename = DataFilePath0 + '/results_2exp_uncontrained_py' + str(100*alpha)+ '.xlsx'    
            workbook = xlsxwriter.Workbook(xlsfilename) # Workbook() takes one, non-optional, argument which is the filename that we want to create.        
            worksheet = workbook.add_worksheet() # The workbook object is then used to add new worksheet via the add_worksheet() method.
            
            # Use the worksheet object to write data via the write() method
            worksheet.write('A1', 'Fname') 
            worksheet.write('B1', 'OBJ')    
            worksheet.write('C1', 'shift') 
            worksheet.write('D1', 'cens')
            worksheet.write('E1', 'lambda1')
            worksheet.write('F1', 'lambda2')
            worksheet.write('G1', 'A1')
            worksheet.write('H1', 'A2')
            worksheet.write('I1', 'K1p')
            worksheet.write('J1', 'K1m')
            worksheet.write('K1', 'K2')
            worksheet.write('L1', 'mRNA')
                
            # for the first line in the table 
            worksheet.write('A2',name) 
            worksheet.write('B2', str(Ominmin))
            worksheet.write(1,2, str(shiftmin)) 
            worksheet.write(1,3, str(cens))
            worksheet.write(1,4, str(round(l1min,4)))
            worksheet.write(1,5, str(round(l2min,4)))
            worksheet.write(1,6, str(round(A1min,4)))
            worksheet.write(1,7, str(round(A2min,4)))
            worksheet.write(1,8, str(round(k1pmin,4)))
            worksheet.write(1,9, str(round(k1mmin,4)))
            worksheet.write(1,10, str(round(k2min,4)))           
            worksheet.write(1,11, 45/mean_w2*60)
    
            # for the second line in the table        
            worksheet.write('A3',name) 
            worksheet.write('B3', '')
            worksheet.write(2,2, '')
            worksheet.write(2,3, '')
            worksheet.write(2,4, str(round(min(l1),4)))
            worksheet.write(2,5, str(round(min(l2),4)))
            worksheet.write(2,6, str(round(min(A1),4)))
            worksheet.write(2,7, str(round(min(A2),4)))
            worksheet.write(2,8, str(round(min(K1p),4)))
            worksheet.write(2,9, str(round(min(K1m),4)))
            worksheet.write(2,10, str(round(min(K2),4)))           

            
            # for the third line in the table
    
            worksheet.write('A3',name) 
            worksheet.write('B3', '')
            worksheet.write(3,2, '')
            worksheet.write(3,3, '')
            worksheet.write(3,4, str(round(max(l1),4)))
            worksheet.write(3,5, str(round(max(l2),4)))
            worksheet.write(3,6, str(round(max(A1),4)))
            worksheet.write(3,7, str(round(max(A2),4)))
            worksheet.write(3,8, str(round(max(K1p),4)))
            worksheet.write(3,9, str(round(max(K1m),4)))
            worksheet.write(3,10, str(round(max(K2),4)))           


        
        workbook.close() # We close the Excel file via the close() method.
    
    
    
    # %%%%%%%%% optimal parameters


    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








