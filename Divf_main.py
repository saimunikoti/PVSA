# -*- coding: utf-8 -*-
"""
Created on Sun Oct 25 16:55:04 2020

@author: saimunikoti
"""
import numpy as np
import scipy.io as io
import matplotlib.pyplot as plt
# from Divf_helping import Divf_simulation

# Ds = Divf_simulation()

#%% initialization
actorlist=[22,14,7,17,8]

phaselist=[2,2,2,1,0,]

## gen covariance matrix of power change

#sigmaS = Ds.get_covariancematrix(actorlist, phaselist, rho_pp=0, rho_qq=0, pvar=4, qvar=1.2, rho_pq=0)

## gene and save power change vector

#Ds.gen_save_powerchangevector(actorlist, phaselist, varp=4, varq=1.2)

#%%  variance of volt change generate deltav_var_actortensor

# deltav_var_actortensor = np.zeros((37,3,len(actorlist)+1)) # 37X3X(# actors +1 )

# ## variance when all actor nodes are present
# fileext = '_'.join(str(elem) for elem in actorlist)

# filename = r"C:\Users\saimunikoti\Manifestation\DOE_Work\PVSA\pvsa_3phase\src\data\Divf\deltav_"+ fileext + "_actor.mat"
# voltchange_tensor= io.loadmat(filename)
# voltchange_tensor = voltchange_tensor['voltchange_tensor']

# deltav_var_actortensor[:,:,0] = np.var(voltchange_tensor, axis=2)

# ## variance when actor nodes var are made zero sequentially
# for countactor in range(0,len(actorlist)):
        
#     tempactorlist = actorlist.copy()
#     tempactorlist.remove(actorlist[countactor])
    
#     fileext = '_'.join(str(elem) for elem in tempactorlist)

#     filename = r"C:\Users\saimunikoti\Manifestation\DOE_Work\PVSA\pvsa_3phase\src\data\Divf\deltav_"+ fileext + "_actor.mat"
#     voltchange_tensor= io.loadmat(filename)
    
#     voltchange_tensor = voltchange_tensor['voltchange_tensor']
        
#     ## get variance of volt change for diff power change
       
#     deltav_var_actortensor[:,:, countactor+1] = np.var(voltchange_tensor, axis=2) # 37X3
    

#%% get difference of variance when var of each actor node is made zero seq
# deltav_vardiff_actortensor = np.zeros((37,3,len(actorlist)))

# for i in range(1,len(actorlist)+1):
#     deltav_vardiff_actortensor[:,:,i-1] = deltav_var_actortensor[:,:,0]- deltav_var_actortensor[:,:,i]
        
#%% find DI nodes in decreasing order

np.where(deltav_vardiff_actortensor[obsnode, phase,:]== np.max(deltav_vardiff_actortensor[obsnode, phase,:]))[0]        

ranked = np.argsort(deltav_vardiff_actortensor, axis=2)
ranked = ranked[:,:,::-1]

actorlistarray = np.array(actorlist)

Diactorstensor = np.zeros((37,3,len(actorlist)))

for countobs in range(37):
    for countph in range(3):
        Diactorstensor[countobs, countph, :] = actorlistarray[ranked[countobs, countph,:]]

for xe, ye in zip(x, y):
    plt.scatter([xe] * len(ye), ye)