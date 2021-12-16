# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 11:34:29 2021

@author: mdouaihy
"""

def constraint0(k):
    constraint =( (sum(k[3:5])>=0) and (k[0]*k[3] + k[1]*k[4] +  k[2]*(1-k[3]-k[4] < 0)) )
    return constraint                                     
