# -*- coding: utf-8 -*-
"""
Created on Sun Oct 25 16:55:04 2020

@author: saimunikoti
"""
import numpy as np
import scipy.io as io
import matplotlib.pyplot as plt
import sys
sys.path.append(".")
from Divf_helping import Divf_simulation, Divf_analytical , Divf_Metrics

Ds = Divf_simulation() 
Da = Divf_analytical(1, 2.7713)
Dm  = Divf_Metrics()

projectpath = r"C:\Users\saimunikoti\Manifestation\DOE_Work\PVSA\pvsa_3phase\src\data\Divf"

#%% initialization
actorlist = [22,14,7,17,8]

phaselist = [2,2,2,1,0,]

varparray = [4,3.5,3.0,2.5,2.0]

varqarray = [1.2,1.1,1.0,0.9,0.8]

rho_pp=0
rho_qq=0
rho_pq=0

################################# Simulation load flwo #############################

#%% gen covariance matrix of power change

sigmaS = Ds.get_covariancematrix(actorlist, phaselist, rho_pp=0, rho_qq=0, pvar=[4,4,4,4,4], qvar=[1.2,1.2,1.2,1.2,1.2], rho_pq=0)

## gene and save power change vector
Ds.gen_save_powerchangevector(actorlist, phaselist, rho_pp, rho_qq, varparray , varqarray, rho_pq )
    
#%% find ranks DI actor nodes with load flow in decreasing order

deltav_var_actortensor, deltav_vardiff_actortensor = Ds.get_deltav_var_actortensor(actorlist, phaselist)
             
DIactorsrank_sim, Diactors_sim = Ds.get_DIranks(deltav_vardiff_actortensor, actorlist)                  

with open(projectpath + "\\DIactorssim_5actors_samevar.pkl", 'wb') as handle:
   
    pickle.dump(Diactors_sim, handle)
    
########################### Theoretical #####################################
#%% ranks and Di actors on the basis of KL distances

voltfilename = projectpath + "\\basevolt_22_14_7_17_8_actor.mat"
basevolt = io.loadmat(voltfilename)
basevolt = basevolt['basevolt'] # kv

sigmaS_Actor = Da.get_sigmaS_A(sigmaS)

DIactorsrank_th, Diactors_th = Da.get_DIactorranks_allobs(basevolt, sigmaS, sigmaS_Actor, actorlist)    

fileext = r"C:\Users\saimunikoti\Manifestation\DOE_Work\PVSA\pvsa_3phase\src\data\Divf"

with open(fileext + "\\DIactorsth_5actors_samevar.pkl", 'wb') as handle:
   
    pickle.dump(Diactors_th, handle)

#%% accuracy
noloadbus = [1,3,4,10,11,15,20,24,25,29,32,33,35]
obsnodearray = np.arange(2,38)
obsnodearray = np.array([elem for elem in obsnodearray if elem not in noloadbus])
 
Acc, TopNacc = Dm.DIaccuracy(Diactors_sim, Diactors_th, obsnodearray, topn, actorlistsize)























    
    


