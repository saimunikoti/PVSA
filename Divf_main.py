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
import pickle

projectpath = r"C:\Users\saimunikoti\Manifestation\DOE_Work\PVSA\pvsa_3phase\src\data\Divf\Actors15_DiffvarSet1\\"

Ds = Divf_simulation(projectpath) 
Da = Divf_analytical(1, 2.7713, projectpath)
Dm  = Divf_Metrics()


#%% initialization
# actorlist = [22,14,7,17,8]
# phaselist = [2, 2, 2,  1,  0]
# varparray = [4, 3.5,3.0, 2.5, 2.0]
# varqarray = [1.2, 1.1, 1.0,  0.9, 0.8]

# actorlist = [7, 8, 14, 17, 18, 22, 27, 28, 31, 36, ]
        
# phaselist = [2, 0, 2,  1,  0,  2,  1,  2,  0 , 1]

# varparray = [4,   4,    3.5, 3.5, 3.0, 3.0, 2.5, 2.5, 2.0, 2.0 ]

# varqarray = [1.2, 1.2,  1.1, 1.1, 1.0, 1.0, 0.9, 0.9, 0.8, 0.8 ]

# 15 actor nodes set 1
actorlist = [7, 8, 9, 12, 14, 17, 18, 22, 26, 27, 28, 30, 31, 34,  36]
        
phaselist = [2, 0, 2,  2,  2,  1,  0,  2,  2,  1,  2,  0,  0,  1  , 1]

varparray = [5, 4.5, 4.0, 3.5, 3.0, 2.5, 2.0, 1.5, 5, 4.5, 4.0, 3.5, 3.0, 2.5, 2.0]

varqarray = [2.2, 2.0, 1.8, 1.6, 1.4, 1.2, 1.0, 2.2, 2.0, 1.8, 1.6, 1.4, 1.2, 1.0, 1.4]

# actor nodes set 2
# actorlist = [7, 8, 9, 12, 14, 17, 19, 22, 26, 27, 28, 31, 34,  36, 37]
        
# phaselist = [2, 0, 2,  2,  2,  1,  0,  2,  2,  1,  2, 0,  1  , 1, 1]

# varparray = [5, 4.5, 4.0, 3.5, 3.0, 2.5, 2.0, 1.5, 5, 4.5, 4.0, 3.5, 3.0, 2.5, 2.0]

# varqarray = [2.2, 2.0, 1.8, 1.6, 1.4, 1.2, 1.0, 2.2, 2.0, 1.8, 1.6, 1.4, 1.2, 1.0, 1.4]


rho_pp=0
rho_qq=0
rho_pq=0

################################# Simulation load flwo #############################

#%% gen covariance matrix of power change

sigmaS = Ds.get_covariancematrix(actorlist, phaselist, rho_pp=0, rho_qq=0, pvar=varparray, qvar=varqarray, rho_pq=0)

## gene and save power change vector
# Ds.gen_save_powerchangevector(actorlist, phaselist, rho_pp, rho_qq, varparray , varqarray, rho_pq )
    
#%% find ranks DI actor nodes with load flow in decreasing order

deltav_var_actortensor, deltav_vardiff_actortensor = Ds.get_deltav_var_actortensor(actorlist, phaselist)
             
DIactorsrank_sim, Diactors_sim = Ds.get_DIranks(deltav_vardiff_actortensor, actorlist)                  

with open(projectpath + "\\DIactorssim_15actors_diffvarset2.pkl", 'wb') as handle:
   
    pickle.dump(Diactors_sim, handle)
    
########################### Theoretical #####################################
#%% ranks and Di actors on the basis of KL distances
fileext = '_'.join(str(elem) for elem in actorlist)

voltfilename = projectpath + "\\basevolt_" + fileext + "_actor.mat"
basevolt = io.loadmat(voltfilename)
basevolt = basevolt['basevolt'] # kv

sigmaS_Actor = Da.get_sigmaS_A(sigmaS)

Diactorsdist_th_kl, Diactors_th_kl , Diactorsdist_th_fd, Diactors_th_fd, \
Diactorsdist_th_bd, Diactors_th_bd, Diactorsdist_th_de, Diactors_th_de = Da.get_DIactorranks_allobs(basevolt, sigmaS, sigmaS_Actor, actorlist)    

#%% saving files
strlist = ['kl','fd','bd','de']

filelist = [Diactors_th_kl, Diactors_th_fd, Diactors_th_bd, Diactors_th_de]

for countfileext, countfile in zip(strlist, filelist):
    
    fileext = "\\DIactorsth_" + str(countfileext)
    
    with open(projectpath + fileext + "_15actors_diffvarset2.pkl", 'wb') as handle:
   
        pickle.dump(countfile, handle)

#%% accuracy and performance
noloadbus = [1,3,4,10,11,15,20,23,24,25,29,32,33,35]
obsnodearray = np.arange(2,38)
obsnodearray = np.array([elem for elem in obsnodearray if elem not in noloadbus])
 
topn = 0.5
actorlistsize= 15
Acc, TopNacc, Di_sim, Di_th_kl, Didist_th_kl = Dm.DIaccuracy(Diactors_sim, Diactorsdist_th_kl, Diactors_th_kl, obsnodearray, topn, actorlist)
Acc, TopNacc, Di_sim, Di_th_fd, Didist_th_fd = Dm.DIaccuracy(Diactors_sim, Diactorsdist_th_fd, Diactors_th_fd, obsnodearray, topn, actorlist)
Acc, TopNacc, Di_sim, Di_th_bd, Didist_th_bd = Dm.DIaccuracy(Diactors_sim, Diactorsdist_th_bd, Diactors_th_bd, obsnodearray, topn, actorlist)
Acc, TopNacc, Di_sim, Di_th_de, Didist_th_de = Dm.DIaccuracy(Diactors_sim, Diactorsdist_th_de,Diactors_th_de, obsnodearray, topn, actorlist)

print(np.mean(TopNacc))

#%% loading saved ranks
with open(projectpath + "\\Actors10_Diffvar\\DIactorssim_10actors_diffvar.pkl", 'rb') as handle:
   
    Diactors_sim = pickle.load(handle)
    
with open(projectpath + "\\Actors10_Diffvar\\DIactorsth_10actors_diffvar.pkl", 'rb') as handle:
   
    Diactors_th = pickle.load(handle)

#%% get obs nodes where rank of actors doesn't match

mismatchnodes_kl = Dm.get_mismatchnodes(Di_sim, Di_th_kl, obsnodearray)
mismatchnodes_fd = Dm.get_mismatchnodes(Di_sim, Di_th_fd, obsnodearray)

mismatchnodes_bd = Dm.get_mismatchnodes(Di_sim, Di_th_bd, obsnodearray)
mismatchnodes_de = Dm.get_mismatchnodes(Di_sim, Di_th_de, obsnodearray)

## voltage influence index for actornode
  
vif_kl, vifnodes_kl = Dm.get_voltinfmetric(Didist_th_kl, actorlist, 3)    
vif_fd, vifnodes_fd = Dm.get_voltinfmetric(Didist_th_fd, actorlist, 3)    
vif_bd, vifnodes_bd = Dm.get_voltinfmetric(Didist_th_bd, actorlist, 3)    

















    
    


