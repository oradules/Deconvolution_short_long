import os
import numpy as np

def fit3_M2(kminmin,Kstoremin,shiftmin,alphamin,censmin,sd):
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

    K3min= -S1;
    K1pmin = 1/2 * ( -L1+S2/S1 + np.sqrt((S1*L1-S2)**2-4*L3*S1)/S1 );
    K2pmin = 1/2 * ( -L1+S2/S1 - np.sqrt((S1*L1-S2)**2-4*L3*S1)/S1 );
    K1mmin = 1/2 * (S1-S2/S1 - (-S1**2*L1+S1*S2+S1*L2-L3+S2**2/S1-S3)/np.sqrt((S1*L1-S2)**2-4*L3*S1));
    K2mmin = 1/2 * (S1-S2/S1 + (-S1**2*L1+S1*S2+S1*L2-L3+S2**2/S1-S3)/np.sqrt((S1*L1-S2)**2-4*L3*S1));
    l1min=l1;l2min=l2;l3min=l3;A1min=A1;A2min=A2;A3min=A3;

    parameters=np.array([K1pmin, K1mmin, K2pmin, K2mmin,K3min])

    indObjAssorted=1
    while np.isnan(parameters).any():
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

        K3min= -S1;        
        K1pmin = 1/2 * ( -L1+S2/S1 + np.sqrt((S1*L1-S2)**2-4*L3*S1)/S1 );
        K2pmin = 1/2 * ( -L1+S2/S1 - np.sqrt((S1*L1-S2)**2-4*L3*S1)/S1 );
        K1mmin = 1/2 * (S1-S2/S1 - (-S1**2*L1+S1*S2+S1*L2-L3+S2**2/S1-S3)/np.sqrt((S1*L1-S2)**2-4*L3*S1));
        K2mmin = 1/2 * (S1-S2/S1 + (-S1**2*L1+S1*S2+S1*L2-L3+S2**2/S1-S3)/np.sqrt((S1*L1-S2)**2-4*L3*S1));
        
        l1min=l1;l2min=l2;l3min=l3;A1min=A1;A2min=A2;A3min=A3;
        indObjAssorted=indObjAssorted+1
        parameters=np.array([K1pmin, K1mmin, K2pmin, K2mmin, K3min])
    
    np.savez('kstore.npz', Kstore=Kstore)

    resM2 = np.array([K1pmin,K1mmin,K2pmin,K2mmin,K3min,l1min,l2min,l3min,A1min,A2min,A3min,Ominmin,45/mean_w2*60,shiftmin,alphamin,censmin,sd[1],sd[0]])

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

    #%%%% model M2
    K3p = -S1;
    K1p = 1/2 * ( -L1+S2/S1 + np.sqrt((S1*L1-S2)**2-4*L3*S1)/S1 );
    K2p = 1/2 * ( -L1+S2/S1 - np.sqrt((S1*L1-S2)**2-4*L3*S1)/S1 );
    K1m = 1/2 * (S1-S2/S1 - (-S1**2*L1+S1*S2+S1*L2-L3+S2**2/S1-S3)/np.sqrt((S1*L1-S2)**2-4*L3*S1));
    K2m = 1/2 * (S1-S2/S1 + (-S1**2*L1+S1*S2+S1*L2-L3+S2**2./S1-S3)/np.sqrt((S1*L1-S2)**2-4*L3*S1));

    KUintervals = np.column_stack((K1p,K1m,K2p,K2m,K3p))
    #validInd = np.where(np.logical_not(np.isnan(KUintervals[:,0:5]).any(axis=1)))
    indobjneg = np.where(np.greater(KUintervals[:,0:5],0).all(axis = 1))[0]
    #######################################################################
    # NS = len(xsgmin) 
    # NL = len(xlmin)

    reslM2= np.array([np.min(K1p[indobjneg]),np.min(K1m[indobjneg]),np.min(K2p[indobjneg]),np.min(K2m[indobjneg]),np.min(K3p[indobjneg]),np.min(l1[indobjneg]),np.min(l2[indobjneg]),np.min(l3[indobjneg]),np.min(A1[indobjneg]),np.min(A2[indobjneg]),np.min(A3[indobjneg])])
    #resl[0:8] = np.max(np.vstack([resl[0:8], np.zeros(8)]), axis=0)
    reshM2= np.array([np.max(K1p[indobjneg]),np.max(K1m[indobjneg]),np.max(K2p[indobjneg]),np.max(K2m[indobjneg]),np.max(K3p[indobjneg]),np.max(l1[indobjneg]),np.max(l2[indobjneg]),np.max(l3[indobjneg]),np.max(A1[indobjneg]),np.max(A2[indobjneg]),np.max(A3[indobjneg])])
    
    return [resM2, reslM2, reshM2]
