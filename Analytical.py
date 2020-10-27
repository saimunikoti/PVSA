# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 22:49:55 2019

@author: saimunikoti
"""

# Calculation of voltage change from the derived expression

import numpy as np
import pandas as pd
import networkx as nx
import pickle
import matplotlib.pyplot as plt
import scipy.io as io
from fitter import Fitter

class Analytical:
    
    def __init__(self, source_node, basekv):
        print("Analytical class is invoked")
        self.source_node = source_node
        self.basekv= basekv
        
    def get_123commonpath(self, actor_node, obs_node) :
        G = nx.Graph()
        edgedata= pd.read_excel(r"C:\Users\saimunikoti\Manifestation\DOE_Work\PVSA\pvsa_3phase\data\raw\feeder123\line data.xls")
        edgedata = edgedata.iloc[2:,:]
        listedges = []
        for countedge in range(edgedata.shape[0]):
            listedges.append((edgedata.iloc[countedge,0], edgedata.iloc[countedge,1]))
        
        G.add_edges_from(listedges)
        G.add_edges_from([(3,1),(13,152),(18,135),(60,160),(61,610),(97,197),(150,149)])
        
        l1 = nx.shortest_path(G, source=150, target= actor_node)
        l2 = nx.shortest_path(G, source=150, target= obs_node)
        
        commonpath = sorted(set(l1) & set(l2), key = l1.index) 
        
        return commonpath

    def get_commonpath(self, actor_node, obs_node):
        
        G = nx.Graph()
#        G.add_edges_from([(799,701),(701,702),(702,703),(702,713),(713,704),
#                          (704,714),(714,718),(704,720),(720,707),(707,724),(707,722),
#                          (720,706),(706,725),(702,705),(705,712),(705,742),(703,727),
#                          (727,744),(744,728),(744,729),(703,730),
#                          (730,709),(709,775),(709,731),(709,708),(708,732),(708,733),(733,734),(734,710),
#                          (710,736),(710,735),(734,737),(737,738),(738,711),(711,740),
#                          (711,741)])
       
        G.add_edges_from([(1,2),(2,3),(3,4),(3,28),(28,29),
                          (29,30),(30,31),(29,32),(32,35),(35,37),(35,36),
                          (32,33),(33,34),(3,25),(25,26),(25,27),(4,5),
                          (5,6),(6,7),(6,8),(4,9),
                          (9,10),(10,24),(10,23),(10,11),(11,12),(11,13),(13,14),(14,15),
                          (15,17),(15,16),(14,18),(18,19),(19,20),(20,21),
                          (20,22)])
        
        l1 = nx.shortest_path(G, source=self.source_node, target= actor_node)
        l2 = nx.shortest_path(G, source=self.source_node, target= obs_node)
        commonpath = sorted(set(l1) & set(l2), key = l1.index) 
        
        return commonpath
    
    def get_impedancematrix(self, y_matrix):
        
        y_bus = y_matrix.iloc[:,1:] # removing bus names from column
        
        for i in range(y_bus.shape[0]):
            
            for col in np.arange(0,y_bus.shape[1],2):
                
                y_bus.iloc[i,col] =complex(y_bus.iloc[i,col], np.double(''.join(y_bus.iloc[i,col+1].split()[1:])))
        
        dropcols=np.arange(1,222,2)
        
        y_bus.drop(y_bus.columns[dropcols],axis=1,inplace=True)
        
        z_bus = np.linalg.inv(y_bus)
        
        z_df = pd.DataFrame(data = z_bus)
        
        bus_number = [int(x) for x in y_matrix.iloc[:,0]]
        
        z_df.insert(0, 'Bus', bus_number)
        
        file = open(r'C:\Users\saimunikoti\Manifestation\DOE_Work\PVSA\pvsa_3phase\data\processed\z_matrix', 'wb')
        # dump information to that file
        pickle.dump(z_df, file)
        file.close()
                    
        
        return z_df
    
    def get_sharedpathimp(self, commonpath):
        
        file = open(r'C:\Users\saimunikoti\Manifestation\DOE_Work\PVSA\pvsa_3phase\data\processed\z_matrix', 'rb')
        # dump information to that file
        z_bus = pickle.load(file)
        # close the file
        file.close()
        
        Total_imp = np.zeros((3,3))
        
        for count in range(len(commonpath)-1):
            
            rowindex = z_bus[z_bus['Bus']==commonpath[count]].index.tolist()
            colindex = z_bus[z_bus['Bus']==commonpath[count+1]].index.tolist()
            
            # plus 1 in colindex because first column in z matrix is bus number
            # plus 3 in row and col index to account all three phases of the buses
            
            path_imp = z_bus.iloc[rowindex[0]:rowindex[0]+3, colindex[0]+1:colindex[0]+3+1] # plus 1 due to bus column in z matrix
            
            Total_imp = (Total_imp) + np.array(path_imp)
               
        return Total_imp
    
    ### shared path impedance for123 netowrk
    
    def get_123pathimpedance(self,actornode, obsnode):
        
        z0 = np.array([[0.00053+1j*0.0,      0.0+1j*0.0,     0.0+1j*0.0], # ohm/mile line 150 to 149
                       [0.0+1j*0.0,       0.00053+1j*0.0,    0.0+1j*0.0],
                       [0.0+1j*0.0,       0.0+1j*0.0,       0.00053 +1j*0.0]]) 
                
    
        z1 = np.array([[ 0.4576+1j*1.0780, 0.1560+1j*0.5017, 0.1535+1j*0.3849],
                       [ 0.1560+1j*0.5017, 0.4666+1j*1.0482, 0.1580+1j*0.4236],
                       [ 0.1535+1j*0.3849, 0.1580+1j*0.4236, 0.4615 +1j*1.0651]])


        z2 = np.array([[ 0.4666+1j*1.0482, 0.1580+1j*0.4236, 0.1560+1j*0.5017],
                       [ 0.1580+1j*0.4236, 0.4615+1j*1.0651, 0.1535+1j*0.3849],
                       [ 0.1560+1j*0.5017,  0.1535+1j*0.3849, 0.4576 +1j*1.0780]])

        
        z3 = np.array([[ 0.4615+1j*1.0651,  0.1535+1j*0.3849, 0.1580+1j*0.4236],
                       [ 0.1535+1j*0.3849, 0.4576+1j*1.0780, 0.1560+1j*0.5017],
                       [ 0.1580+1j*0.4236, 0.1560+1j*0.5017, 0.4666 +1j*1.0482]])
        
        
        z4 = np.array([[ 0.4615+1j*1.0651,  0.1580+1j*0.4236, 0.1535+1j*0.3849],
                       [ 0.1580+1j*0.4236, 0.4666+1j*1.0482, 0.1560+1j*0.5017],
                       [ 0.1535+1j*0.3849, 0.1560+1j*0.5017, 0.4576 +1j*1.0780]])
        
        z5 = np.array([[0.4666+1j*1.0482, 0.1560+1j*0.5017, 0.1580+1j*0.4236],
                       [0.1560+1j*0.5017, 0.4576+1j*1.0780, 0.1535+1j*0.3849],
                       [0.1580+1j*0.4236, 0.1535+1j*0.3849, 0.4615 +1j*1.0651]]) # ohm/mile

        z6 = np.array([[0.4576+1j*1.0780, 0.1535+1j*0.3849, 0.1560+1j*0.5017],
                       [0.1535+1j*0.3849, 0.4615+1j*1.0651, 0.1580+1j*0.4236],
                       [0.1560+1j*0.5017, 0.1580+1j*0.4236, 0.4666 +1j*1.0482]]) # ohm/mile
        
    
        z7 = np.array([[0.4576+1j*1.0780, 0.0+1j*0.0, 0.1535+1j*0.3849],
                       [0.0+1j*0.0,       0.0+1j*0.0, 0.0+1j*0.0],
                       [0.1535+1j*0.3849, 0.0+1j*0.0, 0.4615 +1j*1.0651]]) # ohm/mile


        z8 = np.array([[0.4576+1j*1.0780, 0.1535+1j*0.3849, 0.0+1j*0.0],
                       [0.1535+1j*0.3849, 0.4615+1j*1.0651, 0.0+1j*0.0],
                       [0.0+1j*0.0,       0.0+1j*0.0,       0.0 +1j*0.0]]) # ohm/mile

    
        z9 = np.array([[1.3292+1j*1.3475, 0.0+1j*0.0,       0.0+1j*0.0],
                       [0.0+1j*0.0,       0.0+1j*0.0,       0.0+1j*0.0],
                       [0.0+1j*0.0,       0.0+1j*0.0,       0.0 +1j*0.0]]) 

    
        z10 = np.array([[0.0+1j*0.0,      0.0+1j*0.0,       0.0+1j*0.0],
                       [0.0+1j*0.0,       1.3292+1j*1.3475, 0.0+1j*0.0],
                       [0.0+1j*0.0,       0.0+1j*0.0,       0.0 +1j*0.0]]) 

    
        z11 = np.array([[0.0+1j*0.0,      0.0+1j*0.0,       0.0+1j*0.0],
                       [0.0+1j*0.0,       0.0+1j*0.0,       0.0+1j*0.0],
                       [0.0+1j*0.0,       0.0+1j*0.0,       1.3292 +1j*1.3475]]) 

    
        z12 = np.array([[1.5209+1j*0.7521, 0.5198+1j*0.2775, 0.4924+1j*0.2157],
                       [0.5198+1j*0.2775, 1.5329+1j*1.7162, 0.5198+1j*0.2775],
                       [ 0.4924+1j*0.2157,0.5198+1j*0.2775, 1.5209 +1j*0.7521]])
       
                   
        commonpath = self.get_123commonpath(actornode, obsnode)
        sumimp = 0
        
        linedata = pd.read_excel(r'C:\Users\saimunikoti\Manifestation\DOE_Work\PVSA\pvsa_3phase\data\raw\feeder123\line data.xls') 
        
        for iind in range(len(commonpath)-1):
            
            temprow = linedata[(linedata['Node A']==commonpath[iind]) & (linedata['Node B']==commonpath[iind+1])]
            
            configno = str(temprow.iloc[0,3])
            
            if configno =='0':
               zfactor = z0
               
            elif configno =='1':
               zfactor = z1

            elif configno =='2':
               zfactor = z2

            elif configno =='3':
               zfactor = z3              
               
            elif configno =='4':
               zfactor = z4

            elif configno =='5':
               zfactor = z5  

            elif configno =='6':
               zfactor = z6

            elif configno =='7':
               zfactor = z7                 
               
            elif configno =='8':
               zfactor = z8

            elif configno =='9':
               zfactor = z9

            elif configno =='10':
               zfactor = z10

            elif configno =='11':
               zfactor = z11  

            elif configno =='12':
               zfactor = z12

               
            sumimp = sumimp + (temprow.iloc[0,2]*zfactor)
            
        sumimp = sumimp*0.000189
        
        return sumimp

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
            
    ##### get absolute voltage change at observation node
    def get_voltchange(self, deltap_vec, deltaq_vec, refvolt_actor, actor_node, obs_node):
        
       
        deltas = np.array(([complex(deltap_vec[ind], deltaq_vec[ind]) for ind in range(0,3)] ))
                
#        refvolt_actor = np.array(volt_df.loc[(volt_df['deltap_741']==self.baseload) & (volt_df['Bus']==actor_node), ['pu1','pu2','pu3']])
#        refvolt_actor = refvolt_actor * self.basekv # conversion to volts
                
        Total_imp = self.get_pathimpedance(actor_node, obs_node)
        
        temprightvec = np.transpose(deltas/np.conj(refvolt_actor))
        
        finalvoltage = np.matmul(Total_imp, temprightvec)
        
        finalvoltage = np.absolute(finalvoltage/-2771.3) # pu conversion
        
        finalvoltage = np.reshape(finalvoltage ,(3,)) 
        
        return finalvoltage 
    
    ##### get voltage change at observation node
    
    def get_voltchange_orgsign(self, deltap_vec, deltaq_vec, refvolt_actor, actor_node, obs_node):
    
   
        deltas = np.array(([complex(deltap_vec[ind], deltaq_vec[ind]) for ind in range(0,3)] ))
                            
        Total_imp = self.get_pathimpedance(actor_node, obs_node)
        
        temprightvec = np.transpose(deltas/np.conj(refvolt_actor))
        
        finalvoltage = np.matmul(Total_imp, temprightvec)
        
        finalvoltage = (finalvoltage/-2771.3) # pu conversion
        
        finalvoltage = np.reshape(finalvoltage ,(3,)) 
    
        return finalvoltage 
    
    def plot_deltavoltage(self,simv, analv, phase, ratio, xindices, actor_node, obs_node):
           
           
           xi = list(range(len(xindices))) 
           plt.plot(xi, analv[:,phase], marker='o',color='lightseagreen')
           plt.plot(xi, simv[:,phase], marker='o',color='coral')
           plt.xticks(xi,xindices,fontsize=16)
           plt.yticks(fontsize = 16)
           plt.xlabel('Power change (kw)', fontsize=18)
           plt.ylabel('Magnitude of Voltage chnage (pu)', fontsize=18)
           plt.legend(("Analytical","Load flow"),fontsize=12)
           plt.title("Voltage change at phase " + str(phase+1) +" when actor node is "+str(actor_node)+" and observation node is " + str(obs_node), fontsize=18)
           plt.text(0.9,1,'k= '+str(np.round(ratio[0,phase],2)), transform=plt.gca().transAxes, fontsize=18)
           plt.grid()
           
    def plot_deltavoltsamescale(self, simv, analv, phase, phaseno, xindices, actor_node, obs_node) :
      
           xi = list(range(len(xindices))) 
           plt.scatter(xi, analv[:,phase], marker='o',color='lightseagreen', s=10**2)
           plt.scatter(xi, simv[:,phase], marker='X',color='coral', s=10**2)
           plt.xticks(xi,xindices,fontsize=16)
           plt.yticks(fontsize = 16)
           top = np.max(analv)
           plt.ylim(0, top*1.2)
           SUB = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
           actornode_sub = ('P'+str(actor_node)).translate(SUB)
           obsnode_sub = ('V'+str(obs_node)).translate(SUB)
           plt.xlabel('Power change '+ u'\u0394'+ actornode_sub+' (Kw)', fontsize=18)
           plt.ylabel('Magnitude of Voltage chnage ' + u'\u0394'+ obsnode_sub +' (pu)', fontsize=18)
           plt.legend(("Proposed method","Traditional Load flow"),fontsize=12)
#           plt.title("Voltage change at phase " + str(phase+1) +" when actor node is "+str(actor_node)+" and observation node is " + str(obs_node), fontsize=18)
           plt.grid(True)
           
    ### subplots for voltage change at single observation nodes
                     
    def subplot_deltavolt(temp,simv, analv, xindices, obs_node, phase):
        
        titlename = ['a' ,'b', 'c']
        temp.scatter(xindices, analv[:,phase], marker='o',color='lightseagreen', s=150)
        temp.scatter(xindices, simv[:,phase], marker='X',color='coral', s=150)
        
        plt.xticks(xindices,fontsize=18)
        temp.yaxis.set_tick_params(labelsize=16)
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        top = np.max(analv)
        temp.set_ylim(0, top*1.2)
        temp.set_ylabel('|'+u'\u0394' + 'V| (pu)', fontsize=19)
        temp.legend(("Proposed","Traditional"),fontsize=12,loc=0)
        temp.yaxis.get_offset_text().set_fontsize(16)
        temp.set_title('Phase - '+ titlename[phase], size = 19)
        temp.grid(True)
        
    ## plot delta volt for all nodes

    def subplot_deltavoltallnodes(temp, simv, analv, xindices, phase):
        
        titlename = ['a' ,'b', 'c']
        temp.scatter(xindices, analv[:,phase], marker='o',color='lightseagreen', s=10**2)
        temp.scatter(xindices, simv[:,phase], marker='X',color='coral', s=10**2)
        
        plt.xticks(fontsize=18)
        temp.yaxis.set_tick_params(labelsize=16)
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        top = np.max(analv)
        temp.set_ylim(0, top*1.2)
        temp.set_xlim(1,38)
        temp.set_ylabel('|'+u'\u0394' + 'V| (pu)', fontsize=19)
        temp.legend(("Proposed","Traditional"),fontsize=15)
        temp.yaxis.get_offset_text().set_fontsize(16)
        temp.set_title('Phase - '+ titlename[phase], size = 19)
        temp.grid(True)
        
    def subplotattributes(self, temp, xlabels, ylabels,title):
              
          plt.xticks(fontsize=14)
          temp.yaxis.set_tick_params(labelsize=16)
          temp.xlabel(xlabels, fontsize=16)
          temp.ylabel(ylabels, fontsize=16)
          temp.legend(loc='best', fontsize=12)
          temp.yaxis.get_offset_text().set_fontsize(16)
          temp.set_title(title, size = 19)
          temp.grid(True)
                       
    def check_linearity(self, a, base):
           
           slope = (a[base]-a[base-1])
           c= a[base]-slope*base
           y =[((slope*xind) + c) for xind in range(0,len(a))]
           plt.plot(a,marker='o', color='lightseagreen')
           plt.plot(y,color='coral')
           
           xindices =np.arange(-84, 196, 7) 
           xi = list(range(len(xindices))) 
           plt.plot(xi, a, marker='o',color='lightseagreen')
           plt.plot(xi, y, color='coral')
           plt.xticks(xi,xindices,fontsize=12)
           plt.yticks(fontsize = 12)
           plt.xlabel('Power change (kw)', fontsize=16)
           plt.ylabel('Magnitude of Voltage (pu)', fontsize=16)
           plt.legend(("Load flow value","Linear line"))
           plt.grid()
           
    def check_slope(self, data, actor_node):
           
              t = np.diff(data, axis=0) 
              
              if (t[1,0]>0) or (t[1,1]>0) or (t[1,2]>0):
                  return actor_node
              else:
                     return 'Not found'
                 


class Analytical_randomness():
    
    def __init__(self, networkdim, n_samples, sourcenode, basekv):
        print("class Spatio_temporal_random is invoked")
        self.networkdim = networkdim
        self.n_samples= n_samples
        from Analytical import Analytical
        self.analoutput = Analytical(sourcenode, basekv)
    
    ############## constant shared path impedance matrix- real part    
    def get_CRfunction(self, z, actorbasevolt ):
        temp = (z.real)*np.cos(np.angle(actorbasevolt, deg=False)) - (z.imag)*np.sin(np.angle(actorbasevolt, deg=False))
        cr_p = temp/(np.absolute(actorbasevolt))
        cr_p = cr_p*(-1)
        
        temp = (z.real)*np.sin(np.angle(actorbasevolt, deg=False)) + (z.imag)*np.cos(np.angle(actorbasevolt, deg=False))
        cr_q = temp/(np.absolute(actorbasevolt))
    
        cr = np.concatenate((cr_p, cr_q), axis=0)
    
        return cr
    
    ############## constant shared path impedance matrix- imag part 
    
    def get_CIfunction(self, z, actorbasevolt ):
        temp = (z.real)*np.sin(np.angle(actorbasevolt, deg=False)) + (z.imag)*np.cos(np.angle(actorbasevolt, deg=False))
        ci_p = temp/(np.absolute(actorbasevolt))
        ci_p = ci_p*(-1)
        
        temp = (z.real)*np.cos(np.angle(actorbasevolt, deg=False)) - (z.imag)*np.sin(np.angle(actorbasevolt, deg=False))
        ci_q = temp/(np.absolute(actorbasevolt))
        ci_q = ci_q*(-1)
        
        ci = np.concatenate((ci_p, ci_q), axis=0)
    
        return ci
    
    def get_S(self, countphase, actorlist, actorphase, pvalue, qvalue):
    
        sa=np.zeros(shape=(2*self.networkdim))
        
        for countactor, actor in enumerate(actorlist):
            pindex = actor-1
            qindex = pindex + self.networkdim
            
            if actorphase[countactor]==countphase:
                sa[pindex] = pvalue[countactor]
                sa[qindex] = qvalue[countactor]
        
        return sa 
    
    ########### covariance matrix old
    
    def get_covariancematrix_old(self, actorlist, phaselist, pcov, qcov, pvar, qvar, pqcov):
        
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
                    sigma_s[pindex,qindex] = pqcov
                    sigma_s[qindex,pindex] = pqcov
                    
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
                        sigma_s[pindex,qindex] = pcov
                        sigma_s[qindex,pindex] = pcov
                        
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
                        sigma_s[pindex,qindex] = qcov
                        sigma_s[qindex,pindex] = qcov
                    
        return sigma_s
    ############## covariance matrix of power change       
    def get_covariancematrix(self, actorlist, phaselist, rho_pp, rho_qq, pvar, qvar, rho_pq):
        
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
    
    ##### fit distribution on voltage change
    def fit_dist(self, data):
        f = Fitter(data ,distributions=['nakagami'])
        f.fit()
        nu, loc, omega = f.fitted_param['nakagami']
        return nu,loc,omega
    

    ####### mean of shared path impedance from ieee 37 data for real part of voltage change
    
    def get_Zrstats(self, ztensor):
    
        Zrstat = np.zeros((6*self.networkdim,3), dtype=float)
    
        def get_zrvec(ztensor, phase): 
            part1= (-1)*(ztensor[phase,0].real)
            part1 = np.repeat(part1,37) # last dim for all actor nodes
            
            part2 = ztensor[phase,0].imag
            part2 = np.repeat(part2,37)
            
            part3 = (0.5)*ztensor[phase,1].real - (np.sqrt(3)/2)*ztensor[phase,1].imag
            part3 = np.repeat(part3,37)
            
            part4 = (-1)*(np.sqrt(3)/2)*ztensor[phase,1].real - (1/2)*ztensor[phase,1].imag
            part4 = np.repeat(part4,37)
            
            part5 = (0.5)*ztensor[phase,2].real + (np.sqrt(3)/2)*ztensor[phase,2].imag
            part5 = np.repeat(part5,37)
            
            part6 = (np.sqrt(3)/2)*ztensor[phase,2].real - (1/2)*ztensor[phase,2].imag
            part6 = np.repeat(part6,37)
            
            zrvec = np.concatenate((part1, part2, part3, part4, part5, part6), axis=0)
            
            return zrvec
    
        for countz in range(3): # running for diff phases   
            Zrstat[:,0] = get_zrvec(ztensor, 0) # a phase
            Zrstat[:,1] = get_zrvec(ztensor, 1)
            Zrstat[:,2] = get_zrvec(ztensor, 2)
                   
        return Zrstat 
    
    ####### mean of shared path impedance from ieee 37 data for imag part of voltage change
    def get_Zistats(self, ztensor):
    
        Zistat = np.zeros((6*self.networkdim, 3), dtype=float)
    
        def get_zivec(ztensor, phase): 
            
            part1 = (-1)*ztensor[phase,0].imag # last dim for all actor nodes
            part1 = np.repeat(part1,37)
            
            part2 = (-1)*ztensor[phase,0].real
            part2 = np.repeat(part2,37)
            
            part3 = (np.sqrt(3)/2)*ztensor[phase,1].real + (1/2)*ztensor[phase,1].imag
            part3 = np.repeat(part3,37)
            
            part4 = (1/2)*ztensor[phase,1].real - (np.sqrt(3)/2)*ztensor[phase,1].imag
            part4 = np.repeat(part4,37)
            
            part5 = (-1)*(np.sqrt(3)/2)*ztensor[phase,2].real + (1/2)*ztensor[phase,2].imag
            part5 = np.repeat(part5,37)
            
            part6 = (1/2)*ztensor[phase,2].real + (np.sqrt(3)/2)*ztensor[phase,2].imag
            part6 = np.repeat(part6,37)
            
            zivec = np.concatenate((part1, part2, part3, part4, part5, part6), axis=0)
                        
            return zivec
    
        for countz in range(3): # running for diff phases   
            Zistat[:,0] = get_zivec(ztensor, 0) # a phase
            Zistat[:,1] = get_zivec(ztensor, 1)
            Zistat[:,2] = get_zivec(ztensor, 2)
                   
        return Zistat 
    
    ###### mean and variance of the impedance of the shared path for all actor nodes 
    
    def get_zbase(self):
        ztensor = np.zeros((3,3,37), dtype = complex)

        for actornode in range(37):
            try:
                ztensor[:,:,actornode] = self.analoutput.get_pathimpedance(actornode+1, obsnode=9)
            except:
                ztensor[:,:,actornode] = np.zeros((3,3))
        
        zrmean = np.mean(ztensor.real, axis=2)  
        zxmean = np.mean(ztensor.imag, axis=2)
            
        zrvar = np.var(ztensor.real, axis=2)    
        zxvar = np.var(ztensor.imag, axis=2)   
        
        return zrmean, zxmean, zrvar, zxvar 
    
    ###### covariance of shared path iimpedance matrix
    
    def get_covmatrix_Zr(self, rho_rr, rho_rx, voltphase):
        
        dim = 6*self.networkdim
        
        sigma_z = np.zeros((dim,dim))
        
#        power_index = {}
#        power_index[0]=0
#        power_index[1]=2
#        power_index[2]=4
#        
#        reactivep_index = {}
#        reactivep_index[0]=1
#        reactivep_index[1]=3
#        reactivep_index[2]=5
       
        ######## upadting variance of Zr vector 
        zrmean, zxmean, zrvar, zxvar  = self.get_zbase()
        
#        zvartensor = self.get_Zrstats(zvar)
#        
#        row,col = np.diag_indices(sigma_z.shape[0])
#        sigma_z[row,col] = zvartensor[:,deltavoltphase]
        
        ###### udating covariance 
        
        ### covarinace part1 - resistance part      
        def get_sigmaz1(zrvar,zxvar,sigma_z):
            for i in range(0,37):
                for count,j in enumerate(range(i,222,37)):
                    if count==0:
                        sigma_z[i,j] = zrvar[voltphase,0]
                        
                    elif count== 1:
                        sigma_z[i,j] = (-1)*rho_rx*np.sqrt(zrvar[voltphase,0]*zxvar[voltphase,0])
                        sigma_z[j,i] = (-1)*rho_rx*np.sqrt(zrvar[voltphase,0]*zxvar[voltphase,0])
                        
                    elif count==2:
                        sigma_z[i,j] = (-0.5)*rho_rr*np.sqrt(zrvar[voltphase,0]*zrvar[voltphase,1]) - (-0.866)*rho_rx*np.sqrt(zrvar[voltphase,0]*zxvar[voltphase,1])
                        sigma_z[j,i] = (-0.5)*rho_rr*np.sqrt(zrvar[voltphase,0]*zrvar[voltphase,1]) - (-0.866)*rho_rx*np.sqrt(zrvar[voltphase,0]*zxvar[voltphase,1])
                    
                    elif count ==3:
                        sigma_z[i,j] = (0.866)*rho_rr*np.sqrt(zrvar[voltphase,0]*zrvar[voltphase,1]) -(-0.5)*rho_rx*np.sqrt(zrvar[voltphase,0]*zxvar[voltphase,1])
                        sigma_z[j,i] = (0.866)*rho_rr*np.sqrt(zrvar[voltphase,0]*zrvar[voltphase,1]) -(-0.5)*rho_rx*np.sqrt(zrvar[voltphase,0]*zxvar[voltphase,1])
    
                    elif count==4:
                        sigma_z[i,j] = (-0.5)*rho_rr*np.sqrt(zrvar[voltphase,0]*zrvar[voltphase,2]) + (-0.866)*rho_rx*np.sqrt(zrvar[voltphase,0]*zxvar[voltphase,2])
                        sigma_z[j,i] = (-0.5)*rho_rr*np.sqrt(zrvar[voltphase,0]*zrvar[voltphase,2]) + (-0.866)*rho_rx*np.sqrt(zrvar[voltphase,0]*zxvar[voltphase,2])
    
                    elif count==5:
                        sigma_z[i,j] = (-0.866)*rho_rr*np.sqrt(zrvar[voltphase,0]*zrvar[voltphase,2]) - (-0.5)*rho_rx*np.sqrt(zrvar[voltphase,0]*zxvar[voltphase,2])
                        sigma_z[j,i] = (-0.866)*rho_rr*np.sqrt(zrvar[voltphase,0]*zrvar[voltphase,2]) - (-0.5)*rho_rx*np.sqrt(zrvar[voltphase,0]*zxvar[voltphase,2])
              
            return sigma_z
        
        ### covarinace part - 2 reactance part
        def get_sigmaz2(zrvar,zxvar,sigma_z):
            for i in range(0,37):
                for count,j in enumerate(range(i,222,37)):
                    if count==0:
                        sigma_z[37+i,j] = (-1)*rho_rx*np.sqrt(zxvar[voltphase,0]*zrvar[voltphase,0])
                        sigma_z[j,37+i] = (-1)*rho_rx*np.sqrt(zxvar[voltphase,0]*zrvar[voltphase,0])
                        
                    elif count== 1:
                        sigma_z[37+i,j] = zxvar[voltphase,0]
                        sigma_z[j,37+i] = zxvar[voltphase,0]
                        
                    elif count==2:
                        sigma_z[37+i,j] = (0.5)*rho_rx*np.sqrt(zxvar[voltphase,0]*zrvar[voltphase,1]) - (0.866)*rho_rx*np.sqrt(zxvar[voltphase,0]*zxvar[voltphase,1])
                        sigma_z[j,37+i] = (0.5)*rho_rx*np.sqrt(zxvar[voltphase,0]*zrvar[voltphase,1]) - (0.866)*rho_rx*np.sqrt(zxvar[voltphase,0]*zxvar[voltphase,1])
                    
                    elif count ==3:
                        sigma_z[37+i,j] = (-0.866)*rho_rr*np.sqrt(zxvar[voltphase,0]*zrvar[voltphase,1]) -(0.5)*rho_rx*np.sqrt(zxvar[voltphase,0]*zxvar[voltphase,1])
                        sigma_z[i,37+i] = (-0.866)*rho_rr*np.sqrt(zxvar[voltphase,0]*zrvar[voltphase,1]) -(0.5)*rho_rx*np.sqrt(zxvar[voltphase,0]*zxvar[voltphase,1])
    
                    elif count==4:
                        sigma_z[37+i,j] = (0.5)*rho_rr*np.sqrt(zxvar[voltphase,0]*zrvar[voltphase,2]) + (0.866)*rho_rx*np.sqrt(zxvar[voltphase,0]*zxvar[voltphase,2])
                        sigma_z[i,37+i] = (0.5)*rho_rr*np.sqrt(zxvar[voltphase,0]*zrvar[voltphase,2]) + (0.866)*rho_rx*np.sqrt(zxvar[voltphase,0]*zxvar[voltphase,2])
    
                    elif count==5:
                        sigma_z[37+i,j] = (0.866)*rho_rr*np.sqrt(zxvar[voltphase,0]*zrvar[voltphase,2]) - (0.5)*rho_rx*np.sqrt(zxvar[voltphase,0]*zxvar[voltphase,2])
                        sigma_z[i,37+i] = (0.866)*rho_rr*np.sqrt(zxvar[voltphase,0]*zrvar[voltphase,2]) - (0.5)*rho_rx*np.sqrt(zxvar[voltphase,0]*zxvar[voltphase,2])
                                            
            return sigma_z

        ### covarinace part - 3 resistance and reactance part
        def get_sigmaz3(zrvar,zxvar,sigma_z):        
            for i in range(0,37):
                for count,j in enumerate(range(i,222,37)):
                    if count==0:
                        sigma_z[74+i,j] = (0.5)*rho_rr*np.sqrt(zrvar[voltphase,1]*zrvar[voltphase,0]) - (0.866)*rho_rx*np.sqrt(zxvar[voltphase,1]*zrvar[voltphase,0])
                        sigma_z[j,74+i] = (0.5)*rho_rr*np.sqrt(zrvar[voltphase,1]*zrvar[voltphase,0]) - (0.866)*rho_rx*np.sqrt(zxvar[voltphase,1]*zrvar[voltphase,0])
                        
                    elif count== 1:
                        sigma_z[74+i,j] = (0.5)*rho_rx*np.sqrt(zrvar[voltphase,1]*zxvar[voltphase,0]) - (0.866)*rho_rr*np.sqrt(zxvar[voltphase,1]*zxvar[voltphase,0])
                        sigma_z[j,74+i] = (0.5)*rho_rx*np.sqrt(zrvar[voltphase,1]*zxvar[voltphase,0]) - (0.866)*rho_rr*np.sqrt(zxvar[voltphase,1]*zxvar[voltphase,0])
                                                
                    elif count==2:                        
                        sigma_z[74+i,j] = (0.25)*rho_rr*np.sqrt(zrvar[voltphase,1]*zrvar[voltphase,1]) +(0.25)*rho_rr*np.sqrt(zxvar[voltphase,1]*zxvar[voltphase,1])-(0.866)*rho_rx*np.sqrt(zrvar[voltphase,1]*zxvar[voltphase,1])
                                          
                        sigma_z[j,74+i] = (0.25)*rho_rr*np.sqrt(zrvar[voltphase,1]*zrvar[voltphase,1]) +(0.25)*rho_rr*np.sqrt(zxvar[voltphase,1]*zxvar[voltphase,1]) -(0.866)*rho_rx*np.sqrt(zrvar[voltphase,1]*zxvar[voltphase,1])
                                                        
                    elif count ==3:
                       sigma_z[74+i,j] = (-0.433)*rho_rr*np.sqrt(zrvar[voltphase,1]*zrvar[voltphase,1]) +(0.5)*rho_rx*np.sqrt(zrvar[voltphase,1]*zxvar[voltphase,1])\
                                       +(0.866)*rho_rr*np.sqrt(zxvar[voltphase,1]*zxvar[voltphase,1])
                       sigma_z[j,74+i] = (-0.433)*rho_rr*np.sqrt(zrvar[voltphase,1]*zrvar[voltphase,1]) +(0.5)*rho_rx*np.sqrt(zrvar[voltphase,1]*zxvar[voltphase,1])\
                                       +(0.866)*rho_rr*np.sqrt(zxvar[voltphase,1]*zxvar[voltphase,1])
                   
                    elif count==4:
                       sigma_z[74+i,j] = (0.25)*rho_rr*np.sqrt(zrvar[voltphase,1]*zrvar[voltphase,2]) +(0.433)*rho_rx*np.sqrt(zrvar[voltphase,1]*zxvar[voltphase,2])\
                                      -(0.433)*rho_rx*np.sqrt(zrvar[voltphase,2]*zxvar[voltphase,1]) -(0.75)*rho_rr*np.sqrt(zxvar[voltphase,1]*zxvar[voltphase,2])
                       sigma_z[j,74+i] = (0.25)*rho_rr*np.sqrt(zrvar[voltphase,1]*zrvar[voltphase,2]) +(0.433)*rho_rx*np.sqrt(zrvar[voltphase,1]*zxvar[voltphase,2])\
                                      -(0.433)*rho_rx*np.sqrt(zrvar[voltphase,2]*zxvar[voltphase,1]) -(0.75)*rho_rr*np.sqrt(zxvar[voltphase,1]*zxvar[voltphase,2])
                      
                    elif count==5:
                       sigma_z[74+i,j] = (0.433)*rho_rr*np.sqrt(zrvar[voltphase,1]*zrvar[voltphase,2]) -(0.25)*rho_rx*np.sqrt(zrvar[voltphase,1]*zxvar[voltphase,2])\
                                      -(0.75)*rho_rx*np.sqrt(zrvar[voltphase,2]*zxvar[voltphase,1]) +(0.433)*rho_rr*np.sqrt(zxvar[voltphase,1]*zxvar[voltphase,2])
                       sigma_z[j,74+i] = (0.433)*rho_rr*np.sqrt(zrvar[voltphase,1]*zrvar[voltphase,2]) -(0.25)*rho_rx*np.sqrt(zrvar[voltphase,1]*zxvar[voltphase,2])\
                                      -(0.75)*rho_rx*np.sqrt(zrvar[voltphase,2]*zxvar[voltphase,1]) +(0.433)*rho_rr*np.sqrt(zxvar[voltphase,1]*zxvar[voltphase,2])
                                            
            return sigma_z

        ### covarinace part - 4 resistance and reactance part
        
        def get_sigmaz4(zrvar,zxvar,sigma_z):        
            for i in range(0,37):
                for count,j in enumerate(range(i,222,37)): 
                    if count==0:
                        sigma_z[111+i,j] = (-0.866)*rho_rr*np.sqrt(zrvar[voltphase,1]*zrvar[voltphase,0]) - (0.5)*rho_rx*np.sqrt(zxvar[voltphase,1]*zrvar[voltphase,0])
                        sigma_z[j,111+i] = (-0.866)*rho_rr*np.sqrt(zrvar[voltphase,1]*zrvar[voltphase,0]) - (0.5)*rho_rx*np.sqrt(zxvar[voltphase,1]*zrvar[voltphase,0])
                        
                    elif count== 1:
                        sigma_z[111+i,j] = (-0.866)*rho_rx*np.sqrt(zrvar[voltphase,1]*zxvar[voltphase,0]) - (0.5)*rho_rr*np.sqrt(zxvar[voltphase,1]*zxvar[voltphase,0])
                        sigma_z[j,111+i] = (-0.866)*rho_rx*np.sqrt(zrvar[voltphase,1]*zxvar[voltphase,0]) - (0.5)*rho_rr*np.sqrt(zxvar[voltphase,1]*zxvar[voltphase,0])
                                                
                    elif count==2:
                       sigma_z[111+i,j] = (-0.433)*rho_rr*np.sqrt(zrvar[voltphase,1]*zrvar[voltphase,1]) +(0.5)*rho_rx*np.sqrt(zrvar[voltphase,1]*zxvar[voltphase,1])\
                                       +(0.866)*rho_rr*np.sqrt(zxvar[voltphase,1]*zxvar[voltphase,1])
                       sigma_z[j,111+i] = (-0.433)*rho_rr*np.sqrt(zrvar[voltphase,1]*zrvar[voltphase,1]) +(0.5)*rho_rx*np.sqrt(zrvar[voltphase,1]*zxvar[voltphase,1])\
                                       +(0.866)*rho_rr*np.sqrt(zxvar[voltphase,1]*zxvar[voltphase,1])
                    
                    elif count ==3:
                        sigma_z[111+i,j] = (0.75)*rho_rr*np.sqrt(zrvar[voltphase,1]*zrvar[voltphase,1]) +(0.25)*rho_rr*np.sqrt(zxvar[voltphase,1]*zxvar[voltphase,1])\
                                       +(0.866)*rho_rx*np.sqrt(zrvar[voltphase,1]*zxvar[voltphase,1])
                        sigma_z[j,111+i] = (0.75)*rho_rr*np.sqrt(zrvar[voltphase,1]*zrvar[voltphase,1]) +(0.25)*rho_rr*np.sqrt(zxvar[voltphase,1]*zxvar[voltphase,1])\
                                       +(0.866)*rho_rx*np.sqrt(zrvar[voltphase,1]*zxvar[voltphase,1])
                    
                    elif count==4:
                       sigma_z[111+i,j] = (-0.433)*rho_rr*np.sqrt(zrvar[voltphase,1]*zrvar[voltphase,2]) -(0.75)*rho_rx*np.sqrt(zrvar[voltphase,1]*zxvar[voltphase,2])\
                                      -(0.25)*rho_rx*np.sqrt(zrvar[voltphase,2]*zxvar[voltphase,1]) -(0.433)*rho_rr*np.sqrt(zxvar[voltphase,1]*zxvar[voltphase,2])
                       sigma_z[j,111+i] = (-0.433)*rho_rr*np.sqrt(zrvar[voltphase,1]*zrvar[voltphase,2]) -(0.75)*rho_rx*np.sqrt(zrvar[voltphase,1]*zxvar[voltphase,2])\
                                      -(0.25)*rho_rx*np.sqrt(zrvar[voltphase,2]*zxvar[voltphase,1]) -(0.433)*rho_rr*np.sqrt(zxvar[voltphase,1]*zxvar[voltphase,2])
    
                    elif count==5:
                       sigma_z[111+i,j] = (-0.75)*rho_rr*np.sqrt(zrvar[voltphase,1]*zrvar[voltphase,2]) +(0.433)*rho_rx*np.sqrt(zrvar[voltphase,1]*zxvar[voltphase,2])\
                                      -(0.433)*rho_rx*np.sqrt(zrvar[voltphase,2]*zxvar[voltphase,1]) +(0.25)*rho_rr*np.sqrt(zxvar[voltphase,1]*zxvar[voltphase,2])
                       sigma_z[j,111+i] = (-0.75)*rho_rr*np.sqrt(zrvar[voltphase,1]*zrvar[voltphase,2]) +(0.433)*rho_rx*np.sqrt(zrvar[voltphase,1]*zxvar[voltphase,2])\
                                      -(0.433)*rho_rx*np.sqrt(zrvar[voltphase,2]*zxvar[voltphase,1]) +(0.25)*rho_rr*np.sqrt(zxvar[voltphase,1]*zxvar[voltphase,2])
                                            
            return sigma_z
        
        ### covarinace part - 5 resistance and reactance part
        
        def get_sigmaz5(zrvar,zxvar,sigma_z):        
            for i in range(0,37):
                for count,j in enumerate(range(i,222,37)):
                    if count==0:
                        sigma_z[148+i,j] = (0.5)*rho_rr*np.sqrt(zrvar[voltphase,2]*zrvar[voltphase,0]) + (0.866)*rho_rx*np.sqrt(zxvar[voltphase,2]*zrvar[voltphase,0])
                        sigma_z[j,148+i] = (0.5)*rho_rr*np.sqrt(zrvar[voltphase,2]*zrvar[voltphase,0]) + (0.866)*rho_rx*np.sqrt(zxvar[voltphase,2]*zrvar[voltphase,0])
                        
                    elif count== 1:
                        sigma_z[148+i,j] = (0.5)*rho_rx*np.sqrt(zrvar[voltphase,2]*zxvar[voltphase,0]) + (0.866)*rho_rr*np.sqrt(zxvar[voltphase,2]*zxvar[voltphase,0])
                        sigma_z[j,148+i] = (0.5)*rho_rx*np.sqrt(zrvar[voltphase,2]*zxvar[voltphase,0]) + (0.866)*rho_rr*np.sqrt(zxvar[voltphase,2]*zxvar[voltphase,0])
                                                
                    elif count==2:
                        sigma_z[148+i,j] = (0.25)*rho_rr*np.sqrt(zrvar[voltphase,1]*zrvar[voltphase,2]) +(0.433)*rho_rx*np.sqrt(zrvar[voltphase,1]*zxvar[voltphase,2])\
                                      -(0.433)*rho_rx*np.sqrt(zrvar[voltphase,2]*zxvar[voltphase,1]) -(0.75)*rho_rr*np.sqrt(zxvar[voltphase,1]*zxvar[voltphase,2])
                        sigma_z[j,148+i] = (0.25)*rho_rr*np.sqrt(zrvar[voltphase,1]*zrvar[voltphase,2]) +(0.433)*rho_rx*np.sqrt(zrvar[voltphase,1]*zxvar[voltphase,2])\
                                      -(0.433)*rho_rx*np.sqrt(zrvar[voltphase,2]*zxvar[voltphase,1]) -(0.75)*rho_rr*np.sqrt(zxvar[voltphase,1]*zxvar[voltphase,2])
                    
                    elif count ==3:
                        sigma_z[148+i,j] = (-0.433)*rho_rr*np.sqrt(zrvar[voltphase,1]*zrvar[voltphase,2]) -(0.75)*rho_rx*np.sqrt(zrvar[voltphase,1]*zxvar[voltphase,2])\
                                      -(0.25)*rho_rx*np.sqrt(zrvar[voltphase,2]*zxvar[voltphase,1]) -(0.433)*rho_rr*np.sqrt(zxvar[voltphase,1]*zxvar[voltphase,2])
                        sigma_z[j,148+i] = (-0.433)*rho_rr*np.sqrt(zrvar[voltphase,1]*zrvar[voltphase,2]) -(0.75)*rho_rx*np.sqrt(zrvar[voltphase,1]*zxvar[voltphase,2])\
                                      -(0.25)*rho_rx*np.sqrt(zrvar[voltphase,2]*zxvar[voltphase,1]) -(0.433)*rho_rr*np.sqrt(zxvar[voltphase,1]*zxvar[voltphase,2])
                    
                    elif count==4:
                        sigma_z[148+i,j] = (0.25)*rho_rr*np.sqrt(zrvar[voltphase,2]*zrvar[voltphase,2]) +(0.75)*rho_rr*np.sqrt(zxvar[voltphase,2]*zxvar[voltphase,2])\
                                       +(0.866)*rho_rx*np.sqrt(zrvar[voltphase,2]*zxvar[voltphase,2])
                        sigma_z[j,148+i] = (0.25)*rho_rr*np.sqrt(zrvar[voltphase,2]*zrvar[voltphase,2]) +(0.75)*rho_rr*np.sqrt(zxvar[voltphase,2]*zxvar[voltphase,2])\
                                       +(0.866)*rho_rx*np.sqrt(zrvar[voltphase,2]*zxvar[voltphase,2])
   
                    elif count==5:
                        sigma_z[148+i,j] = (0.433)*rho_rr*np.sqrt(zrvar[voltphase,2]*zrvar[voltphase,2]) -(0.25)*rho_rx*np.sqrt(zrvar[voltphase,2]*zxvar[voltphase,2])\
                                      +(0.75)*rho_rx*np.sqrt(zrvar[voltphase,2]*zxvar[voltphase,2]) -(0.433)*rho_rr*np.sqrt(zxvar[voltphase,2]*zxvar[voltphase,2])
                        sigma_z[j,148+i] = (0.433)*rho_rr*np.sqrt(zrvar[voltphase,2]*zrvar[voltphase,2]) -(0.25)*rho_rx*np.sqrt(zrvar[voltphase,2]*zxvar[voltphase,2])\
                                      +(0.75)*rho_rx*np.sqrt(zrvar[voltphase,2]*zxvar[voltphase,2]) -(0.433)*rho_rr*np.sqrt(zxvar[voltphase,2]*zxvar[voltphase,2])
                                            
            return sigma_z    

        ### covarinace part - 6 resistance and reactance part
        
        def get_sigmaz6(zrvar,zxvar,sigma_z):        
            for i in range(0,37):
                for count,j in enumerate(range(i,222,37)):
                    if count==0:
                        sigma_z[185+i,j] = (-0.866)*rho_rr*np.sqrt(zrvar[voltphase,0]*zrvar[voltphase,2]) - (-0.5)*rho_rx*np.sqrt(zrvar[voltphase,0]*zxvar[voltphase,2])
                        sigma_z[j,185+i] = (-0.866)*rho_rr*np.sqrt(zrvar[voltphase,0]*zrvar[voltphase,2]) - (-0.5)*rho_rx*np.sqrt(zrvar[voltphase,0]*zxvar[voltphase,2])
                        
                    elif count== 1:
                        sigma_z[185+i,j] = (0.866)*rho_rr*np.sqrt(zxvar[voltphase,0]*zrvar[voltphase,2]) - (0.5)*rho_rx*np.sqrt(zxvar[voltphase,0]*zxvar[voltphase,2])
                        sigma_z[j,185+i] = (0.866)*rho_rr*np.sqrt(zxvar[voltphase,0]*zrvar[voltphase,2]) - (0.5)*rho_rx*np.sqrt(zxvar[voltphase,0]*zxvar[voltphase,2])
                                                
                    elif count==2:
                       sigma_z[185+i,j] = (0.433)*rho_rr*np.sqrt(zrvar[voltphase,1]*zrvar[voltphase,2]) -(0.25)*rho_rx*np.sqrt(zrvar[voltphase,1]*zxvar[voltphase,2])\
                                      -(0.75)*rho_rx*np.sqrt(zrvar[voltphase,2]*zxvar[voltphase,1]) +(0.433)*rho_rr*np.sqrt(zxvar[voltphase,1]*zxvar[voltphase,2])
                       sigma_z[j,185+i] = (0.433)*rho_rr*np.sqrt(zrvar[voltphase,1]*zrvar[voltphase,2]) -(0.25)*rho_rx*np.sqrt(zrvar[voltphase,1]*zxvar[voltphase,2])\
                                      -(0.75)*rho_rx*np.sqrt(zrvar[voltphase,2]*zxvar[voltphase,1]) +(0.433)*rho_rr*np.sqrt(zxvar[voltphase,1]*zxvar[voltphase,2])
                    
                    elif count ==3:
                       sigma_z[185+i,j] = (-0.75)*rho_rr*np.sqrt(zrvar[voltphase,1]*zrvar[voltphase,2]) +(0.433)*rho_rx*np.sqrt(zrvar[voltphase,1]*zxvar[voltphase,2])\
                                      -(0.433)*rho_rx*np.sqrt(zrvar[voltphase,2]*zxvar[voltphase,1]) +(0.25)*rho_rr*np.sqrt(zxvar[voltphase,1]*zxvar[voltphase,2])
                       sigma_z[j,185+i] = (-0.75)*rho_rr*np.sqrt(zrvar[voltphase,1]*zrvar[voltphase,2]) +(0.433)*rho_rx*np.sqrt(zrvar[voltphase,1]*zxvar[voltphase,2])\
                                      -(0.433)*rho_rx*np.sqrt(zrvar[voltphase,2]*zxvar[voltphase,1]) +(0.25)*rho_rr*np.sqrt(zxvar[voltphase,1]*zxvar[voltphase,2])
                    
                    elif count==4:
                       sigma_z[185+i,j] = (0.433)*rho_rr*np.sqrt(zrvar[voltphase,2]*zrvar[voltphase,2]) -(0.25)*rho_rx*np.sqrt(zrvar[voltphase,2]*zxvar[voltphase,2])\
                                      +(0.75)*rho_rx*np.sqrt(zrvar[voltphase,2]*zxvar[voltphase,2]) -(0.433)*rho_rr*np.sqrt(zxvar[voltphase,2]*zxvar[voltphase,2])
                       sigma_z[j,185+i] = (0.433)*rho_rr*np.sqrt(zrvar[voltphase,2]*zrvar[voltphase,2]) -(0.25)*rho_rx*np.sqrt(zrvar[voltphase,2]*zxvar[voltphase,2])\
                                      +(0.75)*rho_rx*np.sqrt(zrvar[voltphase,2]*zxvar[voltphase,2]) -(0.433)*rho_rr*np.sqrt(zxvar[voltphase,2]*zxvar[voltphase,2])
    
                    elif count==5:
                        sigma_z[185+i,j] = (0.75)*rho_rr*np.sqrt(zrvar[voltphase,2]*zrvar[voltphase,2]) +(0.25)*rho_rr*np.sqrt(zxvar[voltphase,2]*zxvar[voltphase,2])\
                                       -(0.866)*rho_rx*np.sqrt(zrvar[voltphase,2]*zxvar[voltphase,2])
                        sigma_z[j,185+i] = (0.75)*rho_rr*np.sqrt(zrvar[voltphase,2]*zrvar[voltphase,2]) +(0.25)*rho_rr*np.sqrt(zxvar[voltphase,2]*zxvar[voltphase,2])\
                                       -(0.866)*rho_rx*np.sqrt(zrvar[voltphase,2]*zxvar[voltphase,2])
                                            
            return sigma_z 
        
        sigma_z = np.zeros((dim,dim))
        
        sigma_1= get_sigmaz1(zrvar,zxvar,sigma_z)  
        sigma_2= get_sigmaz2(zrvar,zxvar,sigma_z)   
        sigma_3= get_sigmaz3(zrvar,zxvar,sigma_z)   
        sigma_4= get_sigmaz4(zrvar,zxvar,sigma_z)   
        sigma_5= get_sigmaz5(zrvar,zxvar,sigma_z)   
        sigma_6= get_sigmaz6(zrvar,zxvar,sigma_z)  
        
        sigma_z = sigma_1+sigma_2+sigma_3+sigma_4+sigma_5+sigma_6
         
        return sigma_z

    ###### random vector of shared path impedance - real part
    
    def get_Zr(self, ztensor_array):
    
        Zr = np.zeros((6*self.networkdim, self.n_samples,3), dtype=float)
    
        def get_zrvec(ztensor, phase):   
            part1 = (-1)*ztensor[phase,0,:].real # last dim for all actor nodes
            part2 = ztensor[phase,0,:].imag
            
            part3 = (0.5)*ztensor[phase,1,:].real - (np.sqrt(3)/2)*ztensor[phase,1,:].imag
            part4 = (-1)*(np.sqrt(3)/2)*ztensor[phase,1,:].real - (1/2)*ztensor[phase,1,:].imag
            
            part5 = (0.5)*ztensor[phase,2,:].real + (np.sqrt(3)/2)*ztensor[phase,2,:].imag
            part6 = (np.sqrt(3)/2)*ztensor[phase,2,:].real - (1/2)*ztensor[phase,2,:].imag
            
            zrvec = np.concatenate((part1, part2, part3, part4, part5, part6), axis=0)
            
            return zrvec
    
        for countz, ztensor in enumerate(ztensor_array): # running for diff cases   
            Zr[:,countz,0] = get_zrvec(ztensor, 0) # diff phases
            Zr[:,countz,1] = get_zrvec(ztensor, 1)
            Zr[:,countz,2]= get_zrvec(ztensor, 2)
                   
        return Zr
    
    ###### random vector of shared path impedance matrix- imag part
    
    def get_Zi(self, ztensor_array):
        
        Zi = np.zeros((6*self.networkdim, self.n_samples,3), dtype=float)
        
        def get_zivec(ztensor, phase):   
            part1 = (-1)*ztensor[phase,0,:].imag # last dim for all actor nodes
            part2 = (-1)*ztensor[phase,0,:].real
            
            part3 = (np.sqrt(3)/2)*ztensor[phase,1,:].real + (1/2)*ztensor[phase,1,:].imag
            part4 = (1/2)*ztensor[phase,1,:].real - (np.sqrt(3)/2)*ztensor[phase,1,:].imag
            
            part5 = (-1)*(np.sqrt(3)/2)*ztensor[phase,2,:].real + (1/2)*ztensor[phase,2,:].imag
            part6 = (1/2)*ztensor[phase,2,:].real + (np.sqrt(3)/2)*ztensor[phase,2,:].imag
            
            zivec = np.concatenate((part1, part2, part3, part4, part5, part6), axis=0)
            
            return zivec
        
        for countz, ztensor in enumerate(ztensor_array): # running for diff cases   
            Zi[:,countz,0] = get_zivec(ztensor, 0) # diff phases
            Zi[:,countz,1] = get_zivec(ztensor, 1)
            Zi[:,countz,2]= get_zivec(ztensor, 2)
                   
        return Zi
    
    ###### cobvariance matrix for delta power change 
    
    def get_covmatrix_spatiotemp(self, actorlist, phaselist, pcov, qcov, pvar, qvar, pqcov, basevolt):
        
        sigma_s = self.get_covariancematrix(actorlist, phaselist, pcov, qcov, pvar, qvar, pqcov)
        
        temp_basevolt = np.concatenate((basevolt, basevolt), axis=0)    

        temp_volt = np.concatenate((temp_basevolt[:,0], temp_basevolt[:,1], temp_basevolt[:,2]), axis=0)        
        
        sigmaS = sigma_s/temp_volt
        
        return sigmaS
                     
        
        
        
        
        
                  
           
           
           
           
           
           
           
           
           

           
        


    
