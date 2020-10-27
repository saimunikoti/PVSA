# -*- coding: utf-8 -*-
"""
Analtical and MC simulation helping functions

"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import scipy.io as io

class Divf_analytical():
    
    
    def __init__(self):
        print("Divf analytical class is invoked")
        self.actornodes = np.arange(1, 37)
        self.networkdim =37
        

    def get_CRfunction(self, z, actorbasevolt ):
        
        """
        input z = Raa or Rab or Rac (NX1 or 1X1)
        
        input basevolt: (NX1 or 1X1)
        
        output cr: (2NX1 or 2X1)
            
        """
        
        temp = (z.real)*np.cos(np.angle(actorbasevolt, deg=False)) - (z.imag)*np.sin(np.angle(actorbasevolt, deg=False))
        cr_p = temp/(np.absolute(actorbasevolt))
        cr_p = cr_p*(-1)
        
        temp = (z.real)*np.sin(np.angle(actorbasevolt, deg=False)) + (z.imag)*np.cos(np.angle(actorbasevolt, deg=False))
        cr_q = temp/(np.absolute(actorbasevolt))
        cr_q = cr_q*(-1)
        
        cr = np.concatenate((cr_p, cr_q), axis=0)

        return cr 
    

    def get_CIfunction(self, z, actorbasevolt ):
            
        """
        input z = Raa or Rab or Rac (NX1 or 1X1)
        
        input basevolt: (NX1 or 1X1)
        
        output ci: (2NX1 or 2X1)
            
        """
        temp = (z.real)*np.sin(np.angle(actorbasevolt, deg=False)) + (z.imag)*np.cos(np.angle(actorbasevolt, deg=False))
        ci_p = temp/(np.absolute(actorbasevolt))
        ci_p = ci_p*(-1)
        
        temp = (z.imag)*np.sin(np.angle(actorbasevolt, deg=False)) - (z.real)*np.cos(np.angle(actorbasevolt, deg=False))
        ci_q = temp/(np.absolute(actorbasevolt))
        ci_q = ci_q*(-1)
        
        ci = np.concatenate((ci_p, ci_q), axis=0)
    
        return ci  
    
    def compute_CR(self, ztensor, basevolt, voltphase):
        
        """
        input voltphase  = a/b/c phase for which volt change is calculated 
        
        output cr_a for real part : (6NX1)
        
        """

        cr_a = np.array([])
        for countphase in range(3):
            cr_temp = self.get_CRfunction(ztensor[voltphase,countphase,:], basevolt[:,countphase] )
            cr_a = np.concatenate((cr_a, cr_temp), axis=0)
            
        cr_a = np.reshape(cr_a,(cr_a.shape[0],1))   # returns in kilo 
        return cr_a    
    
    
    def compute_CI(self,ztensor, basevolt , voltphase):
        """
        input voltphase  = a/b/c phase for which volt change is calculated 
        
        output ci_a for imag part : 6NX1
        
        """
        ci_a = np.array([])
        for countphase in range(3):
            
            ci_temp = self.get_CIfunction(ztensor[voltphase,countphase,:], basevolt[:,countphase] )
            ci_a = np.concatenate((ci_a, ci_temp), axis=0)
            
        ci_a = np.reshape(ci_a,(ci_a.shape[0],1))   # returns in kilo    
        return ci_a      
    
    ############ get impedance of the shared path between obs node and actor. 
    
    def get_pathimpedance(self,actornode, obsnode):

        z721 = np.array([ [0.2926+1j*0.1973,   0.0673-1j*0.0368,   0.0337-1j*0.0417],
         [0.0673-1j*0.0368 ,  0.2646+1j*0.1900 ,  0.0673-1j*0.0368],
         [0.0337-1j*0.0417  , 0.0673-1j*0.0368  , 0.2926+1j*0.1973] ])# ohm per mile
          
        z722 = np.array([ [0.4751+1j*0.2973 ,  0.1629-1j*0.0326,   0.1234-1j*0.0607],
         [0.1629-1j*0.0326   ,0.4488+1j*0.2678 ,  0.1629-1j*0.0326],
         [0.1234-1j*0.0607 ,  0.1629-1j*0.0326 ,  0.4751+1j*0.2973 ]])# ohm per mile
     
        z723 = np.array([[ 1.2936+1j*0.6713 ,  0.4871+1j*0.2111,   0.4585+1j*0.1521],
         [0.4871+1j*0.2111  , 1.3022+1j*0.6326   ,0.4871+1j*0.2111],
         [0.4585+1j*0.1521 ,  0.4871+1j*0.2111,   1.2936+1j*0.6713] ])# ohm per mile
     
        z724 = np.array([ [2.0952+1j*0.7758 ,  0.5204+1j*0.2738 ,  0.4926+1j*0.2123],
         [0.5204+1j*0.2738  , 2.1068+1j*0.7398  , 0.5204+1j*0.2738],
         [0.4926+1j*0.2123 ,  0.5204+1j*0.2738 ,  2.0952+1j*0.7758] ])# ohm per mile
       
               
        commonpath = self.get_commonpath(actornode, obsnode)
        sumimp = 0
        
        linedata = pd.read_excel(r'C:\Users\saimunikoti\Manifestation\DOE_Work\PVSA\pvsa_3phase\data\processed\Lineconfig.xlsx') 
        for iind in range(len(commonpath)-1):
            
            temprow = linedata[(linedata['Node A']==commonpath[iind]) & (linedata['Node B']==commonpath[iind+1])]
            
            configno = str(temprow.iloc[0,3])
            
            if configno =='721':
               zfactor = z721
               
            elif configno =='722':
               zfactor = z722

            elif configno =='723':
               zfactor = z723

            elif configno =='724':
               zfactor = z724               
               
            sumimp = sumimp + (temprow.iloc[0,2]*zfactor)
            
        sumimp = sumimp*0.000189
        
        return sumimp
    
    def get_ztensor(self, obsnode):
        """
        input obsnode: the node number in original format
        output- ztensor: 3X3X37
        
        """
        
        ztensor = np.zeros((3,3,self.networkdim), dtype = complex)
       
        for actornode in range(self.networkdim):
            
            try:
                ztensor[:,:,actornode] = self.get_pathimpedance(actornode+1, obsnode)
            except:
                ztensor[:,:,actornode] = np.zeros((3,3))
                
        return ztensor
    
    # get covariance matrix of power change vector for each single actor node
    def get_sigmaS_A(self, sigmaS):
        
        """
        input - sigmaS: 6Nx6N
        
        output- sigmaS_actor: 6X6X37
        
        """
        # covariance matrix of an actor node for all the three phases
        sigmaS_Actor = np.zeros((6,6, self.networkdim))
        
        for actornode in range(1, self.networkdim+1):   
            
            ## fill diag elements
        
            for i in range(6):            
                sigmaS_Actor[i,i, actornode-1] = sigmaS[self.networkdim*i + actornode-1, self.networkdim*i + actornode-1]
                                                                              
            ## fill non-diag eleents
            sigmaS_Actor[0,1, actornode-1] = sigmaS[self.networkdim*0 + actornode-1, self.networkdim*1 + actornode-1]
            sigmaS_Actor[1,0, actornode-1] =  sigmaS_Actor[0,1, actornode-1]
        
            sigmaS_Actor[2,3, actornode-1] = sigmaS[self.networkdim*2 + actornode -1 , self.networkdim*3 + actornode-1]
            sigmaS_Actor[3,2, actornode-1] =  sigmaS_Actor[2,3, actornode-1]
        
            sigmaS_Actor[4,5, actornode-1] = sigmaS[self.networkdim*4 + actornode -1,  self.networkdim*5 + actornode -1]
            sigmaS_Actor[5,4, actornode-1] =  sigmaS_Actor[4,5, actornode-1]
              
        return sigmaS_Actor  
    
    # probability distribution of voolt change due to all actor nodes                        
    def gen_voltchngdist_multipleactor(self, ztensor, sigmaS, basevolt, nsamples):
        """
        Input : ztensor: 3x3x37
                sigmaS : 6Nx6N
                basevolt: 37X3
            
        Output: sigma_deltavolt_allphase: 2X2X3 covariance matrix of real and 
                imaginary part of voltage change due to all actor
        """
        dvoltreal_abc= np.zeros((nsamples,3))
        dvoltimag_abc= np.zeros((nsamples,3))
        
        # covariance matrix for all the three phases
        sigma_deltavolt_allphase = np.zeros((2,2,3))
        CR_allphase = np.zeros((6*self.networkdim,3))
        CI_allphase = np.zeros((6*self.networkdim,3))
        
        for voltphase in range(3): # loop for volt chnage in three phases

            CR = self.compute_CR(ztensor, basevolt, voltphase) # 6NX1
            CI = self.compute_CI(ztensor, basevolt, voltphase) # 6NX1
            
            CR_allphase[:,voltphase] = CR
            CI_allphase[:,voltphase] = CI
            
            tempm = np.matmul(sigmaS,CR)  
            sigma_deltavolt_allphase[0,0,voltphase] = np.matmul(CR.T, tempm)
            dvoltreal_abc[:,voltphase] = np.random.multivariate_normal([0],sigma_deltavolt_allphase[0,0,voltphase], nsamples)[:,0]  
    
            tempm = np.matmul(sigmaS,CI)  
            sigma_deltavolt_allphase[1,1, voltphase] = np.matmul(CI.T, tempm)
            dvoltimag_abc[:,voltphase] = np.random.multivariate_normal([0], sigma_deltavolt_allphase[1,1, voltphase], nsamples)[:,0]   

            tempm = np.matmul(sigmaS,CI)  
            sigma_deltavolt_allphase[0,1, voltphase] = np.matmul(CR.T, tempm)
            sigma_deltavolt_allphase[1,0, voltphase] = sigma_deltavolt_allphase[0,1, voltphase]  
            
              
        dvoltreal_abc = dvoltreal_abc/2771.3  # volt pu
        dvoltimag_abc = dvoltimag_abc/2771.3  # volt
    
        dvoltabs_abc = np.sqrt((dvoltreal_abc)**2 + (dvoltimag_abc)**2)

        return sigma_deltavolt_allphase, CR_allphase, CI_allphase, dvoltreal_abc, dvoltimag_abc , dvoltabs_abc 


    ## probability dist. of volt change due to single actor node        
    def gen_voltchngdist_singleactor(self, CR_allphase, CI_allphase, sigmaS_Actor, actornode, nsamples=100000):
        
        """
        input - CR_allphase: 6Nx3
                CI_allphase: 6Nx3
                sigmaS_Actor: 6X6X37
                actornode: actornode number out of 37
                nsamples : samples for voltage change distribution
                
        output- sigma_deltavolt_Actor_allphase: 2X2X3 cov matrix of voltage 
            change gaussian random vector due to single actor
        
        """

        dvoltreal_abc= np.zeros((nsamples,3))
        dvoltimag_abc= np.zeros((nsamples,3))
        sigma_deltavolt_Actor_allphase = np.zeros((2,2,3))
        
        # loop for dist. of volt chnage in three phases                    
        for voltphase in range(3): 

            CR_A = [CR_allphase[self.networkdim*i + actornode-1, voltphase] for i in range(6)]
            CI_A = [CI_allphase[self.networkdim*i + actornode-1, voltphase] for i in range(6)]
                        
            tempm = np.matmul(sigmaS_Actor[:,:,actornode-1],CR_A)  
            sigma_deltavolt_Actor_allphase[0,0,voltphase] = np.matmul(CR_A.T, tempm)
            dvoltreal_abc[:,voltphase] = np.random.multivariate_normal([0], sigma_deltavolt_Actor_allphase[0,0,voltphase], nsamples)[:,0]  
    
            tempm = np.matmul(sigmaS_Actor[:,:,actornode-1],CI_A)  
            sigma_deltavolt_Actor_allphase[1,1, voltphase] = np.matmul(CI_A.T, tempm)
            dvoltimag_abc[:,voltphase] = np.random.multivariate_normal([0], sigma_deltavolt_Actor_allphase[1,1, voltphase], nsamples)[:,0]   

            tempm = np.matmul(sigmaS_Actor[:,:,actornode-1],CI_A)  
            sigma_deltavolt_Actor_allphase[0,1, voltphase] = np.matmul(CR_A.T, tempm)
            sigma_deltavolt_Actor_allphase[1,0, voltphase] = sigma_deltavolt_Actor_allphase[0,1, voltphase] 
            
            
        dvoltreal_abc = dvoltreal_abc/2771.3  # volt pu
        dvoltimag_abc = dvoltimag_abc/2771.3  # volt
    
        dvoltabs_abc = np.sqrt((dvoltreal_abc)**2 + (dvoltimag_abc)**2)

        return sigma_deltavolt_Actor_allphase, dvoltreal_abc, dvoltimag_abc, dvoltabs_abc
    
    def gen_voltchngdist_singleactor_tensor(self, CR_allphase, CI_allphase, sigmaS_Actor, nsamples=100000):
        """
        input - CR_allphase, CI_allphase, sigmaS_Actor, # of actornodes
                
        output- collect volt change dustribution individually for each actor
        """
        sigma_deltavolt_Actor_allphase_tensor = np.zeros((self.networkdim,2,2,3))
        
        for actornode in range(1, self.networkdim+1):
            sigma_deltavolt_Actor_allphase_tensor[actornode-1,:,:,:] = self.gen_voltchngdist_singleactor(CR_allphase, CI_allphase, sigmaS_Actor, actornode)
        
        return sigma_deltavolt_Actor_allphase_tensor
    
    ## get optimal actor nodes for KL divergence metric
    def get_optimnode_kl(self, obsnode, sigma_deltavolt_allphase, sigma_deltavolt_Actor_allphase_tensor):
        """
        input -sigma_deltavolt_allphase: 2X2X3
            sigma_deltavolt_Actor_allphase_tensor: 37X2X2X3
                
        output- DI_KL: vector of distances for all actor and phases (37,3)
        """
        
        Di_KL = np.zeros(shape=(self.networkdim,3))
        
        ## zero mean power change 
        mu = np.zeros(shape=(2,1))
        mu_A = np.zeros((2,1))
        
        for actornode in range(1, self.networkdim+1): 
            
            for voltphase in range(3):
                
                sigma_dvoltinv = np.inv(sigma_deltavolt_allphase[:,:,voltphase])
                sigma_dvoltA = sigma_deltavolt_Actor_allphase_tensor[actornode-1,:,:,voltphase]
            
                Di_KL[actornode-1,voltphase] = 0.5(np.trace(np.matmul(sigma_dvoltinv, sigma_dvoltA)) -2 + \
                               np.matmul( (mu-mu_A), np.matmul( sigma_dvoltinv, (mu-mu_A)) ) + \
                               np.log(np.linalg.det(sigma_deltavolt_allphase[:,:, voltphase])/np.linalg.det(sigma_dvoltA)) )
            
        return Di_KL # check the correct numbering of actor nodes
            
    def get_optimnode_fd(self, obsnode, phase):
        
        sigma_dvolt = self.gen_voltchngdist_multipleactor(obsnode, phase)
        
        DI_fd = np.zeros(shape=(len(self.actornodes),1))
        for countactor in self.actornodes:
            sigma_dvoltA = self.gen_voltchngdist_singleactor(obsnode, countactor, phase)
            
            DI_fd[countactor] = "Frechet formula"
            
        return np.argmin(DI_fd) # check the correct numbering of actor nodes

##########################  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##################
## class with help functions for generating distribution from Load flow simulation

class Divf_simulation():
        
    def __init__(self):
        
        print("Divf simulation class is invoked")
        
        self.networkdim =37
    
    ## generate covariance matrixc for given variance and correlation coefficient
    def get_covariancematrix(self, actorlist, phaselist, rho_pp, rho_qq, pvar, qvar, rho_pq):
        
        """
        input- actorlist and corresponding phaselist, variances and 
               correlation coefficient
               
        output- sigmaS matrix (6NX6N)
        """
        dim = 6*self.networkdim
        
        sigma_s = np.zeros((dim,dim))
        power_index = {}
        power_index[0]=0
        power_index[1]=2
        power_index[2]=4
        
        reactivep_index = {}
        reactivep_index[0]=1
        reactivep_index[1]=3
        reactivep_index[2]=5
               
        ######### updating cross cov between p and q 
        for (activep,phasep) in zip(actorlist, phaselist):
            phasenop = power_index[phasep]
            pindex= phasenop*self.networkdim + (activep-1)
            for (reactivep,phaseq) in zip(actorlist, phaselist):

                if phasep == phaseq:
                    phasenoq = reactivep_index[phaseq]
                    qindex =  phasenoq*self.networkdim + (reactivep-1)
                    sigma_s[pindex,qindex] = rho_pq*np.sqrt(pvar)*np.sqrt(qvar)
                    sigma_s[qindex,pindex] = rho_pq*np.sqrt(pvar)*np.sqrt(qvar)
                    
        del pindex, qindex, phasenop, phasenoq
           
        ######### updating cov between active power
        for (activep,phasep) in zip(actorlist, phaselist):
            
            ############## active power cov
            phasenop = power_index[phasep]
            pindex = phasenop*self.networkdim+(activep-1)
            for (reactivep,phaseq) in zip(actorlist, phaselist):

                if phasep == phaseq: # active p cov for same phase 
                    phasenoq = power_index[phaseq]
                    qindex =  phasenoq*self.networkdim+(reactivep-1)
                    if activep == reactivep: # active p for variance 
                        sigma_s[pindex,qindex] = pvar
                        sigma_s[qindex,pindex] = pvar
                    else:    # active p cov but same phase 
                        sigma_s[pindex,qindex] = rho_pp*pvar
                        sigma_s[qindex,pindex] = rho_pp*pvar
                        
            del pindex, qindex, phasenop, phasenoq ,reactivep ,phaseq 
            
            ############## reactive power cov
            phasenop = reactivep_index[phasep]
            pindex = phasenop*self.networkdim+(activep-1)            
            for (reactivep,phaseq) in zip(actorlist, phaselist):

                if phasep == phaseq: # reactive q cov for same phase 
                    phasenoq = reactivep_index[phaseq]
                    qindex =  phasenoq*self.networkdim+(reactivep-1)
                    if activep == reactivep: # reactive q for variance 
                        sigma_s[pindex,qindex] = qvar
                        sigma_s[qindex,pindex] = qvar
                    else:    # reactive q cov but same phase 
                        sigma_s[pindex,qindex] = rho_qq*qvar
                        sigma_s[qindex,pindex] = rho_qq*qvar
                    
        return sigma_s 
    
    ## generate and save power change vector         
    def gen_save_powerchangevector(self, actorlist, phaselist, varp=4, varq =1.2, nsamples=100000):
        
        """
        input: sigmaS cov matrix, actorlist, phaselist, varp, varq
        output: multiple instances of power change vector (nsamples X 222) with all actor nodes
                multiple instnces of power change vector with vari of each actor nodes as zero
        
        """ 
        
        sigmaS = self.get_covariancematrix(actorlist, phaselist, rho_pp=0, rho_qq=0, pvar=varp, qvar=varq, rho_pq=0)
        
       
        mean = np.zeros(shape=(6*self.networkdim))        
        powervec = np.random.multivariate_normal(mean, sigmaS, nsamples)
        
        fileext = '_'.join(str(elem) for elem in actorlist)
        tempdic={}
        tempdic['powervec'] = powervec
        filename = r"C:\Users\saimunikoti\Manifestation\DOE_Work\PVSA\pvsa_3phase\src\data\Divf\powervec_"+ fileext + ".mat"
        io.savemat(filename, tempdic)
        
        ## gen and save power change vector with variance set to zero 
        
        for countactor in range(0,len(actorlist)):
            
            tempactorlist = actorlist.copy()
            tempphaselist = phaselist.copy()
            tempactorlist.remove(actorlist[countactor])
            tempphaselist.remove(phaselist[countactor])
                    
            sigmaStemp = self.get_covariancematrix(tempactorlist, tempphaselist, rho_pp=0, rho_qq=0, pvar=4, qvar=1.2, rho_pq=0)
                             
            powervec = np.random.multivariate_normal(mean, sigmaStemp, nsamples)
            
            fileext = '_'.join(str(elem) for elem in tempactorlist)
            
            tempdic={}
            tempdic['powervec'] = powervec
            
            filename = r"C:\Users\saimunikoti\Manifestation\DOE_Work\PVSA\pvsa_3phase\src\data\Divf\powervec_"+ fileext + ".mat"
            
            io.savemat(filename, tempdic)         
    
    def get_deltav_var_actortensor(self,actorlist, phaselist):
        
        """
        input: actorlist, phaselist
        output: variance of all nodes with all phases. Third dimension is the 
                different scenarios of voltage change with variance of each 
                actor as zero. deltav_var_actortensor : (37,3, #actor +1)
                
                deltav_vardiff_actortensor (37,3,#actors): difference of variance with all
                actor nodes to individual cases with each actor node var as 0
        """ 
        
        deltav_var_actortensor = np.zeros((37,3,len(actorlist)+1)) # 37X3X(# actors +1 )

        ## variance when all actor nodes are present
        fileext = '_'.join(str(elem) for elem in actorlist)
        
        filename = r"C:\Users\saimunikoti\Manifestation\DOE_Work\PVSA\pvsa_3phase\src\data\Divf\deltav_"+ fileext + "_actor.mat"
        voltchange_tensor= io.loadmat(filename)
        voltchange_tensor = voltchange_tensor['voltchange_tensor']
        
        deltav_var_actortensor[:,:,0] = np.var(voltchange_tensor, axis=2)
        
        ## variance when actor nodes var are made zero sequentially
        for countactor in range(0,len(actorlist)):
                
            tempactorlist = actorlist.copy()
            tempactorlist.remove(actorlist[countactor])
            
            fileext = '_'.join(str(elem) for elem in tempactorlist)
        
            filename = r"C:\Users\saimunikoti\Manifestation\DOE_Work\PVSA\pvsa_3phase\src\data\Divf\deltav_"+ fileext + "_actor.mat"
            voltchange_tensor= io.loadmat(filename)
            
            voltchange_tensor = voltchange_tensor['voltchange_tensor']
                
            ## get variance of volt change for diff power change
               
            deltav_var_actortensor[:,:, countactor+1] = np.var(voltchange_tensor, axis=2) # 37X3
        
        deltav_vardiff_actortensor = np.zeros((37,3,len(actorlist)))

        for i in range(1,len(actorlist)+1):
            deltav_vardiff_actortensor[:,:,i-1] = deltav_var_actortensor[:,:,0]- deltav_var_actortensor[:,:,i]

        return deltav_var_actortensor, deltav_vardiff_actortensor
        
    ## find DI nodes rank in decreasing order
    def get_DIranks(self, deltav_vardiff_actortensor, actorlist) :
        
        """
        
        input: deltav_vardiff_actortensor: tensor with varianes differences  
        output: ranked tensor with ranks in decreasing order (37,3,#actor)
                                                              
                actors in the order of their ranks for each observation node 
                Diactors_tensor : (37,3,#actors) 
                                              
        """ 
        
        ranked = np.argsort(deltav_vardiff_actortensor, axis=2)
        ranked = ranked[:,:,::-1]
        
        actorlistarray = np.array(actorlist)
        Diactors_tensor = np.zeros((self.networkdim,3,len(actorlist)))

        for countobs in range(self.networkdim):
            for countph in range(3):
                Diactors_tensor[countobs, countph, :] = actorlistarray[ranked[countobs, countph,:]]

        return ranked, Diactors_tensor
       
        
        
        
        
        