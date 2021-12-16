# -*- coding: utf-8 -*-
"""
Created on Fri Jun  4 12:13:27 2021

@author: mdouaihy
"""

import numpy as np
from scipy.io import savemat
from scipy.io import loadmat
import os
import shutil

def npzToMat(path,name,outputpath):#careful name here is without .npz
    content=np.load(path+name+'.npz')
    d = dict(zip(("{}".format(k) for k in content), (content[k] for k in content)))
    fname=name+'.mat'
    savemat(outputpath+fname,d)
    
    
def MatToNpz(path,name,outputpath):#careful name here is without .npz
    content=loadmat(path+name+'.mat')
    
    
    
  #  d = dict(zip(("{}".format(k) for k in content), (content[k] for k in content)))
    fname=name+'.npz'
    savez_dict = dict()
    
    name_list=list(content)[3:]
    for i in name_list:
        savez_dict[i] = content[i] 
    np.savez(outputpath+fname,**savez_dict)
    
def ListingFiles(DataFilePath,extension):
    DataFilePath=DataFilePath
    file_name_list = np.array(os.listdir(DataFilePath)) # list of subdirectories containing data from different genotypes
    nexp = len(file_name_list) # number of files inside he folder

    #just to avoid erros in the opening of the files
    file_name_list_modified=np.array([])
    for i in range(nexp):
        if  extension in file_name_list[i]:
            file_name_list_modified=np.append(file_name_list_modified,file_name_list[i])
                
    file_name_list=file_name_list_modified            
    file_name_list=np.unique(file_name_list)
    nexp=len(file_name_list)
    return file_name_list, nexp
    
    

#pathh='/home/rachel/Documents/Manuscript/deconvolutionResults/resultDec/fitVar.mat'
pathh='/home/rachel/Documents/longMovie/long_movie_artificial/variables/'

outputpath='/home/rachel/Documents/longMovie/long_movie_artificial/variables/'
extension='.mat'
lists,nexp=ListingFiles(pathh,extension)
path=pathh
outputpath=outputpath

# if os.path.exists(outputpath): 
#     shutil.rmtree(outputpath, ignore_errors = True)  
# os.mkdir(outputpath)
# for i in range(len(framelen)):
#     for x in range(1,6):
#         path=pathh+str(framelen[i]) +'/'
#         outputpath=path
#         name=namee+str(framelen[i])+'_'+str(x)
#         npzToMat(path,name,outputpath)

for i in range(nexp):
    name=lists[i].replace(extension,'')
    print(name)
    if extension=='.mat':
        MatToNpz(path, name, outputpath)
    elif extension=='.npz':
        npzToMat(path, name, outputpath)
        