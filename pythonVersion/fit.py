#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 24 00:34:49 2021

@author: rachel
"""

import numpy as np
import os
import shutil
import matplotlib.pyplot as plt
import math
import pandas as pd
from extentForPlot import *
from scipy.io import loadmat
import seaborn as sns
from posPred2cPosPred import *
from ecdfEstimate import ecdf_bounds
import math
import random
from constraint0 import constraint0
from interpolateMetm import interpolateMetm
from fitSumOf2ExponentialsMain import fitSumOf2ExponentialsMain
#from fitSumOf3ExponentialsMain import *
from fitSumOf3ExponentialsContrained import fitSumOf3ExponentialsContrained
from scipy import interpolate
from scipy.optimize import least_squares
from fit2 import *
from fit3 import *
#from common_fit2_part import fit2
#from common_fit3_part import fit3


class fit:
    def __init__(self,path,parameter,combined,visualize,outputpath):
        self.path=path
        self.parameterpath=parameter
        self.visualize = visualize
        ### parameters used
        filecontent=np.load(self.parameterpath)
        Polym_speed=filecontent['Polym_speed']
        TaillePreMarq=filecontent['TaillePreMarq']
        TailleSeqMarq=filecontent['TailleSeqMarq']
        TaillePostMarq=filecontent['TaillePostMarq']
        EspaceInterPolyMin=filecontent['EspaceInterPolyMin']
        FrameLen=filecontent['FrameLen']
        Intensity_for_1_Polym=filecontent['Intensity_for_1_Polym']
        FreqEchImg=filecontent['FreqEchImg']
        DureeSignal=filecontent['DureeSignal']
        
        FreqEchSimu = 1/(EspaceInterPolyMin/Polym_speed) # how many interval(possible poly start position) in 1s
        self.FreqEchSimu = FreqEchSimu 
        
        ####### parameters for the plots
        fsz=16 #figure size
        outliers_short = 1
        outliers_long = 1
         ## function needed to set the the parameters for the color map
        cm_jet= plt.cm.get_cmap('jet') # set the colormap to jet array
        
        DataFilePath0 = outputpath+'Results'
        if os.path.exists(DataFilePath0):
            shutil.rmtree(DataFilePath0, ignore_errors = True)

        os.mkdir(DataFilePath0)


        ### creating of xls sheet for all kind of models
        
        ## 2 states
        
        xlsfilename2states = DataFilePath0 + '/fit2_results.xlsx'
        # OBJ	shift	alpha	cens	lambda1	lambda2	A1	A2	k2	k1m	k1p	mRNA
        # Setting the Names of the Data output in the excel sheet
        xls_cont = pd.DataFrame({'Data': [],'k1p': [], 'k1m': [], 'k2': [],
                                 'l1':[],'l2':[],'A1':[],'A2':[], 'Obj': [],'mRNA':[],'shift': [],
                                 'alpha': [],'cens': [],'samples': [],'Frames': []})    

        # Create a Pandas Excel writer using XlsxWriter as the engine.
        writer2states = pd.ExcelWriter(xlsfilename2states, engine='xlsxwriter')

        # Convert the dataframe to an XlsxWriter Excel object.
        xls_cont.to_excel(writer2states, sheet_name='Sheet1', index=False)
        
        #############################
        
        ### 3 states M1
        xlsfilename3statesM1 = DataFilePath0 + '/fit3M1_results.xlsx'
         
        # Setting the Names of the Data output in the excel sheet
        xls_cont = pd.DataFrame({'Data': [],'k1p' : [],'k1m' : [],'k2p' : [],
                                 'k2m' : [],'k3' : [],'lambda1': [],'lambda2': [],
                                 'lambda3': [],'A1': [],'A2': [],'A3': [],
                                 'Obj': [],'mRNA':[],'shift': [],
                                 'alpha': [],'cens': [],'samples': [],'Frames': []})  
        
        
        # Create a Pandas Excel writer using XlsxWriter as the engine.
        writer3statesM1 = pd.ExcelWriter(xlsfilename3statesM1, engine='xlsxwriter')

        # Convert the dataframe to an XlsxWriter Excel object.
        xls_cont.to_excel(writer3statesM1, sheet_name='Sheet1', index=False)
        
        #########################
        
        ### 3 states M2
        xlsfilename3statesM2 = DataFilePath0 + '/fit3M2_results.xlsx'
         
        # Setting the Names of the Data output in the excel sheet
        xls_cont = pd.DataFrame({'Data': [],'k1p' : [],'k1m' : [],'k2p' : [],
                                 'k2m' : [],'k3' : [],'lambda1': [],'lambda2': [],
                                 'lambda3': [],'A1': [],'A2': [],'A3': [],'S1': [],
                                 'Obj': [],'mRNA':[],'shift': [],
                                 'alpha': [],'cens': [],'samples': [],'Frames': []}) 
        
        
        # Create a Pandas Excel writer using XlsxWriter as the engine.
        writer3statesM2 = pd.ExcelWriter(xlsfilename3statesM2, engine='xlsxwriter')

        # Convert the dataframe to an XlsxWriter Excel object.
        xls_cont.to_excel(writer3statesM2, sheet_name='Sheet1', index=False)
        
        #########################

        
        ####################################################################
        
        ## loading the result of the deconvolution
        NPZFilePath = path;  
        file_name_list = np.array(os.listdir(NPZFilePath)) # list of the data
        # only consider the names without 'long' in the name
        #file_name_list = np.array(os.listdir('./resultDec/')) # list of the data
        longFiles = list(filter(lambda x:'long' in x, file_name_list))
        file_name_list = np.setdiff1d(file_name_list, longFiles, assume_unique=True)
        self.file_name_list=file_name_list
        nexp = len(file_name_list) # length of the list
        nfiles=nexp
        #######################################################################
        pooled =combined; #### if this is 1, pool all the result files from NPZfilePath
        ############### if not, use each file separately########################
    
        if not pooled:
            nfiles=nexp
        else:
            nfiles=1        
            
    
        ### starting the fit for each file
        for ifile in range(nfiles):

            if pooled: 
                
            ### first compute max dimension
                nmax=1
                nmaxpos=1

                for iii in range(nexp):
                    fname=file_name_list[iii]
                    ffname=NPZFilePath+fname
                    if '.npz' in ffname:
                        fnameContent=np.load(ffname)
                    else:
                        fnameContent=loadmat(ffname)
                    DataPred = fnameContent['DataPred']
                    DataExp=fnameContent['DataExp']
                    Fit = fnameContent['Fit']
                    PosPred=fnameContent['PosPred']
                    n2 = DataExp.shape
                    n3 = PosPred.shape
            
                    if n2[0] >nmax:
                        nmax = n2[0]
            
                    if n3[0]> nmaxpos:
                        nmaxpos = n3[0]
               
                #### lump files, Procustes method 
                dataExp=np.empty((nmax, 0), int)
                dataPred=np.empty((nmax,0), int)
                posPred= np.empty((nmaxpos,0), int)
                tmax = np.empty((0, n2[1]), int)

            
                for iii in range(nexp): 
                    fname=file_name_list[iii]
                    ffname=NPZFilePath+fname
                    if '.npz' in ffname:
                        fnameContent=np.load(ffname)
                    else:
                        fnameContent=loadmat(ffname)
                    DataPred = fnameContent['DataPred']
                    DataExp=fnameContent['DataExp']
                    Fit = fnameContent['Fit']
                    PosPred=fnameContent['PosPred']
                    cPosPred = np.asarray(posPred2cPosPred(PosPred))
                    n2 = DataExp.shape
                    n3 = PosPred.shape
      
                    DataExp=np.append(DataExp,np.zeros((nmax-n2[0],n2[1])),axis=0)
                    DataPred=np.append(DataPred,np.zeros((nmax-n2[0],n2[1])),axis=0)
                    PosPred=np.append(PosPred,np.zeros((nmaxpos-n3[0],n3[1])),axis=0)
            
                    # we are adding all the data from different files together
                    dataExp = np.append(dataExp, DataExp, axis=1) 
                    dataPred = np.append(dataPred, DataPred, axis=1)
                    posPred = np.append(posPred, PosPred, axis=1)

                    tmax=np.append(tmax, n2[0]/FreqEchImg*np.ones(n2[1]))
            
        
                DataExp = dataExp.copy()
                DataPred=dataPred.copy()
                PosPred = posPred.copy()
            else:
                fname = file_name_list[ifile]
                print(fname)
                #### full path file name
                ffname = NPZFilePath + fname
                if '.npz' in ffname:
                    fnameContent=np.load(ffname)
                else:
                    fnameContent=loadmat(ffname)
                DataPred = fnameContent['DataPred']
                DataExp=fnameContent['DataExp']
                Fit = fnameContent['Fit']
                PosPred=fnameContent['PosPred']
                cPosPred = np.asarray(posPred2cPosPred(PosPred))
                n2=DataExp.shape
                n_cells = len(DataPred[0])
                nSamples = n2[1]
                frame_num=n2[0] ### number of frames
                DureeSimu = frame_num*FrameLen  ### film duration in s
                DureeAnalysee = DureeSignal + DureeSimu ###(s)
                tmax=DureeAnalysee ### duration of timeseries (consider all the polymerases between 0 and tmax)
#                tmax=n2[0]/FreqEchImg*np.ones(n2[1]) #### movie length, the same for all nuclei in a data sets 
                
                # Load long movie data
                if '.npz' in fname:
                    iend=fname.index('.npz')
                else:
                    iend=fname.index('.mat') 
                lname=fname[0:iend] 
                lname=lname.split('_')[1]
                lfname = list(filter(lambda x:lname in x, longFiles))
                ffname = NPZFilePath + lfname[0]
                if '.npz' in ffname:
                    fnameContent=np.load(ffname)
                else:
                    fnameContent=loadmat(ffname)            
                DataExpLong = fnameContent['DataExpLong']
                DataExpLongRaw = fnameContent['DataExpLongRaw']
                Time = fnameContent['Time']
                Time = Time.reshape(np.size(Time),)
            
            ### extract short name from result file name
            if '.npz' in fname:
                iend=fname.index('.npz')
            else:
                iend=fname.index('.mat') 
            name=fname[0:iend] 
            name=name.split('_')[1]
            self.name= name 
            
            
            ### where to write figure files 
            dirwrite = DataFilePath0+'/'+name+'_result'
            if os.path.exists(dirwrite):
                shutil.rmtree(dirwrite, ignore_errors = True)

            os.mkdir(dirwrite)

            n = DataExp.shape
            nexp = n[1]
            
            ## parameters
            DureeSimu = n[0]*FrameLen #in s
            frame_num = n[0]
            DureeAnalysee = DureeSignal + DureeSimu # (s)
            num_possible_poly = round(DureeAnalysee/(EspaceInterPolyMin/Polym_speed))

                            ### Figure showing Data Signal Prediction
            h=plt.figure(40)
            sz= DataPred.shape
            Y_normal = np.arange(1,sz[1]+1)
            Y=Y_normal[::-1]
            X = np.arange(0, sz[0])/FreqEchImg/60
            plt.imshow(DataPred.T, cmap=cm_jet, extent=extentForPlot(X).result + extentForPlot(Y).result, aspect='auto', origin='upper')
            plt.xlabel('Time [min]', fontsize=12)
            plt.ylabel('Transcription site', fontsize=12)
            cb= plt.colorbar()
            cb.ax.tick_params(labelsize=fsz)
            figfile=dirwrite+'/DataPred_'+name+'.pdf'
            h.savefig(figfile, dpi=800) 
            plt.close()

            ### Figure showing Data Signal Experimental
            h = plt.figure(50)
            plt.imshow(DataExp.T, cmap=cm_jet, extent=extentForPlot(X).result + extentForPlot(Y).result, aspect='auto', origin='upper')
            plt.xlabel('Time [min]', fontsize=12)
            plt.ylabel('Transcription site', fontsize=12)
            plt.colorbar()
            figfile=dirwrite+'/DataExp_'+name+'.pdf'
            h.savefig(figfile, dpi=800) 
            plt.close()

            ### Figure showing Data Position Prediction
            h=plt.figure(60)
            Y_normal=np.arange(1, len(PosPred[0])+1)
            Y=Y_normal[::-1]
            X=np.arange(0,len(PosPred))*EspaceInterPolyMin/Polym_speed/60  -(TaillePreMarq+TailleSeqMarq+TaillePostMarq)/Polym_speed/60 ### time
            plt.imshow(PosPred.T, cmap='gray', extent=extentForPlot(X).result + extentForPlot(Y).result, aspect='auto', origin='upper')
            plt.xlabel('Time [min]', fontsize=12)
            plt.ylabel('Transcription site', fontsize=12)
            figfile=dirwrite+'/PosPred'+name+'.pdf'
            h.savefig(figfile, dpi=800) 
            plt.close()

           
            ######### eliminate outliers handling short movie ########


#             """"Several observed transcription sites had abnormal behaviour (too many or too few events).
#             The decision was made to take them off the data set as follows which is what we are doing with 
#             outliners handling"""" 

            if outliers_short: 

                Ev = np.zeros((n_cells)) ####  based on quartile of the numbers of polymerases which is equivalent to saying: Compute the amount of events of transcription that happened during the movie
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
            ### compute distribution of spacings
            nn=cPosPred.shape

            dt=np.array([])
            dtc=np.array([])
            tmax=0
            
            figfile=dirwrite +'/PosPred'+name+'.txt' ###  text file for pol positions 
            fid = open(figfile,'w+')

            
            for i in range(nn[1]): #for all cells
                timesline = (np.where(PosPred[:,i]==1)[0]+1) / FreqEchSimu -(TaillePreMarq+TailleSeqMarq+TaillePostMarq)/Polym_speed 
                fid.writelines([' \n'+ str(timesline/60)])
                indices = cPosPred[0,i][0] #### indice is the estimate position of polymerase for each cell
                times = indices / FreqEchSimu ### last transcription in the movie
                Mtimes=max(times)

                if Mtimes > tmax:
                    tmax = Mtimes

            fid.close()           

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

            sd=DataExp.shape

            ###### eliminate outliers handling long movie #######
            #         '''
            # """""""""" Some transcription sites in long movies data set also show unusual behaviours being active (FI=0) or inactive
            # (FI=1) during the entire movie. We exclude these outliers as we did it for the short movies, but based on the
            # fraction of inactivity for each transcription site.""""""""""""""""'''
        
            
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
            tstep = 3
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
            lDEL = len(DataExpLong[0])
            
            [res, resl, resh] = fit2(dirwrite,name,dt,dtg,censored,censored_short,wt,wtc,lDEL,Total,Ninactive,visualize,time,sd)
            df1 = pd.DataFrame([res.tolist(), resl.tolist(), resh.tolist()]) 
            df1.to_excel(writer2states,sheet_name='Sheet1', startrow=4*(ifile+1)-3, startcol=1, header=False, index=False)
            df2 = pd.DataFrame([name]) #filename
            df2.to_excel(writer2states,sheet_name='Sheet1', startrow=4*(1+ifile)-3, startcol=0, header=False, index=False)            
        
            ########## save parameters results for 3 state model
            [resM1, reslM1, reshM1]=fit3(dirwrite,name,dt,dtg,censored,censored_short,wt,wtc,lDEL,Total,Ninactive,visualize,time,sd)
            
            #### Model M1            
            df1M1 = pd.DataFrame([resM1.tolist(), #best result
                    reslM1.tolist(), # low 
                    reshM1.tolist()]) # high
            
            df1M1.to_excel(writer3statesM1,sheet_name='Sheet1', startrow=4*(ifile+1)-3, startcol=1, header=False, index=False)
            df2M1 = pd.DataFrame([name.replace('result_','')]) #filename
            df2M1.to_excel(writer3statesM1,sheet_name='Sheet1', startrow=4*(ifile+1)-3, startcol=0, header=False, index=False)

        writer2states.save()
        writer3statesM1.save()