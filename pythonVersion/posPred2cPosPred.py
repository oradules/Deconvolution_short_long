#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 19 14:37:23 2021

@author: rachel
"""
import numpy as np
def posPred2cPosPred(PosPred):
    cPosPred = np.zeros((1,np.shape(PosPred)[1]), dtype=object)
    for ii in range(np.shape(PosPred)[1]):
        posic = (np.where(PosPred[:,ii]==1)[0])
        posic = posic.reshape((1,np.size(posic)))
        cPosPred[0][ii] = posic
    return cPosPred




