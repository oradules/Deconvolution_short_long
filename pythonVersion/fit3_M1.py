import os
import numpy as np

def fit3_M1(kminmin,Kstoremin,shiftmin,alphamin,censmin,sd):
    kmin=kminmin[:];
    Kstore=Kstoremin[:];
    shiftmin;alphamin;censmin;sd[1];sd[0];
    
    ### finding the optimal finite parameters
    ind = np.where(np.max(np.abs(Kstore[:,0:3].imag )  ,axis=1) <1e-10)[0] #takes the part where we have the imaginary part of lambda i almost 0

    ### this will reshape the matrix store in a way that we have an assending order of the last column objective 
    StoreSorted=Kstore[np.argsort(Kstore[ind, -1])]
    objmin = np.min(Kstore[ind,-1]).real #minimum of the objectives for all the previous lambdas

    indmin=np.argmin(Kstore[ind,-1]) #finding the positions of the lambda where we have the lowest objective
    imin = ind[indmin] #finding the index of the positions of the lambda where we have the lowest objective

    kmin = Kstore[imin,0:5].real# taking the lambda i's and A i's that are the fittest
    Ominmin=Kstore[imin,-1]

    l1=kmin[0];
    l2=kmin[1];
    l3=kmin[2];
    A1=kmin[3];
    A2=kmin[4];
    A3=1-A1-A2;
    mean_w2 = -A1/l1-A2/l2-A3/l3; #%%%% average waiting time from parameters
    #mean_w1=mw1opt;
    L1=l1+l2+l3;
    L2=(np.multiply(l1,l2)+np.multiply(l1,l3)+np.multiply(l2,l3))
    L3=l1*l2*l3;
    S1=A1*l1+A2*l2+A3*l3;
    S2=A1*l1**2+A2*l2**2+A3*l3**2;
    S3=A1*l1**3+A2*l2**3+A3*l3**3;

    k1min=-L3*(S1**2-S2)/(S2**2-S1*S3);
    k2min=-(S2**2-S1*S3)/S1/(S1**2-S2);
    k3min=-S1;
    k4min=(S1**2-S2)/S1;
    k5min=-A1*A2*A3*(l1-l2)**2*(l1-l3)**2*(l2-l3)**2*S1/(S1**2-S2)/(S2**2-S1*S3);
    
    l1min=l1;l2min=l2;l3min=l3;A1min=A1;A2min=A2;A3min=A3;

    parameters=np.array([k1min, k5min, k2min, k4min, k3min])

    indObjAssorted=1
    while np.isnan(parameters).any():
        #print(parameters)
        objmin = StoreSorted[indObjAssorted,-1] # np.min(store[ind,-1]).real #minimum of the objectives for all the previous lambdas

        imin = indObjAssorted# ind[indmin] #finding the index of the positions of the lambda where we have the lowest objective

        kmin = StoreSorted[imin,0:5].real# taking the lambda i's and A i's that are the fittest
        Ominmin = StoreSorted[imin,-1]
        # compute 5 rates k1p,m k2p,m k3 from the 5 parameters
        l1=kmin[0];
        l2=kmin[1];
        l3=kmin[2];
        A1=kmin[3];
        A2=kmin[4];
        A3=1-A1-A2;
        mean_w2 = -A1/l1-A2/l2-A3/l3; #%%%% average waiting time from parameters
        #mean_w1=mw1opt;
        L1=l1+l2+l3;
        L2=(np.multiply(l1,l2)+np.multiply(l1,l3)+np.multiply(l2,l3))
        L3=l1*l2*l3;
        S1=A1*l1+A2*l2+A3*l3;
        S2=A1*l1**2+A2*l2**2+A3*l3**2;
        S3=A1*l1**3+A2*l2**3+A3*l3**3;

        k1min=-L3*(S1**2-S2)/(S2**2-S1*S3);
        k2min=-(S2**2-S1*S3)/S1/(S1**2-S2);
        k3min=-S1;
        k4min=(S1**2-S2)/S1;
        k5min=-A1*A2*A3*(l1-l2)**2*(l1-l3)**2*(l2-l3)**2*S1/(S1**2-S2)/(S2**2-S1*S3);
        
        l1min=l1;l2min=l2;l3min=l3;A1min=A1;A2min=A2;A3min=A3;
        indObjAssorted=indObjAssorted+1
        parameters=np.array([k1min, k5min, k2min, k4min, K3min])
    np.savez('kstore.npz', Kstore=Kstore)

    resM1 = np.array([k1min,k5min,k2min,k4min,k3min,l1min,l2min,l3min,A1min,A2min,A3min,Ominmin,45/mean_w2*60,shiftmin,alphamin,censmin,sd[1],sd[0]])

    #%%%%%%%%%%%%%%%% compute parameters with error bars %%%%%%%%%%%%%%
    ###################################################################
    overflow=1;# %%%%% 100% overflow

    #%%%% select near-optimal parameters
    #%%%% O between Ominmin and Ominmin*(1+overflow)  %%%  
    ind = np.where( (Kstore[:,-1] < (1+overflow)*objmin) & (np.max(np.abs(Kstore[:,0:3].imag   ),axis=1) <1e-10)) [0]
    Ksel = Kstore[ind,0:5].real
    #Ksel = Kstore[ Kstore[ind:,5] < Ominmin*(1+overflow) , :  ]
    #%%%% compute intervals for parameters

    l1=Ksel[:,0];
    l2=Ksel[:,1];
    l3=Ksel[:,2];
    A1=Ksel[:,3];
    A2=Ksel[:,4];
    A3=1-A1-A2;

    MRNA = 45/(-A1/l1-A2/l2-A3/l3)*60;

    L1=l1+l2+l3;
    L2=l1*l2+l1*l3+l2*l3;
    L3=l1*l2*l3;
    S1=A1*l1+A2*l2+A3*l3;
    S2=A1*l1**2+A2*l2**2+A3*l3**2;
    S3=A1*l1**3+A2*l2**3+A3*l3**3;

    #%%%% model M1
    K1=-L3*(S1**2-S2)/(S2**2-S1*S3); #%%% k1p
    K2=-(S2**2-S1*S3)/S1/(S1**2-S2); #%%% k2p
    K3=-S1; #%%%% k3
    K4=(S1**2-S2)/S1; #%%% k2m
    K5=-A1*A2*A3*(l1-l2)**2*(l1-l3)**2*(l2-l3)**2*S1/(S1**2-S2)/(S2**2-S1*S3); #%%% k1m

    KUintervals = np.column_stack((K1,K5,K2,K4,K3))
    #validInd = np.where(np.logical_not(np.isnan(KUintervals[:,0:5]).any(axis=1)))
    indobjneg = np.where(np.greater(KUintervals[:,0:5],0).all(axis = 1))[0]
    #######################################################################
    # NS = len(xsgmin) 
    # NL = len(xlmin)

    reslM1= np.array([np.min(K1[indobjneg]),np.min(K5[indobjneg]),np.min(K2[indobjneg]),np.min(K4[indobjneg]),np.min(K3[indobjneg]),np.min(l1[indobjneg]),np.min(l2[indobjneg]),np.min(l3[indobjneg]),np.min(A1[indobjneg]),np.min(A2[indobjneg]),np.min(A3[indobjneg])])
    #resl[0:8] = np.max(np.vstack([resl[0:8], np.zeros(8)]), axis=0)
    reshM1= np.array([np.max(K1[indobjneg]),np.max(K5[indobjneg]),np.max(K2[indobjneg]),np.max(K4[indobjneg]),np.max(K3[indobjneg]),np.max(l1[indobjneg]),np.max(l2[indobjneg]),np.max(l3[indobjneg]),np.max(A1[indobjneg]),np.max(A2[indobjneg]),np.max(A3[indobjneg])])
    
    return [resM1, reslM1, reshM1]
