# -*- coding: utf-8 -*-
"""
Created on Wed Aug  8 09:55:38 2018

@author: smh
"""
import os
from datetime import datetime
from Parameter import ParameterList
from STEPSRiverNetwork import STEPSRiverNetwork
from InputData import InputData
import numpy as np

class RunFactory(object):
    
    def __init__(self,fpath,fname):
        """
        fdir (string): path of project folder
        """
#        self.message("create runfactory: " + fdir)
        self.runlist = ParameterList(fpath,fname,",")
        self.runs = []
        
 

    def setup(self):
        """
        """        
        # create runs
        for rinfo in self.runlist: 
            if rinfo.simulation == True:
                
                # read input data
                inpD = InputData(os.path.join(rinfo.fpath,rinfo.key))
                cmp = inpD.SubstanceList[rinfo.substance][0]
                
                # get reach üarameter as numpy list
                DENS = np.array([i.dens for i in inpD.ReachList])
                POROSITY = np.array([i.porosity for i in inpD.ReachList])
                OC = np.array([i.oc for i in inpD.ReachList])
                DEPTH_SED = np.array([i.depth_sed for i in inpD.ReachList])
                DEPTH_SED_DEEP = np.array([i.depth_sed_deep for i in inpD.ReachList])
                
                # setup STEPSRivernetowrk
                has_par = [rinfo.MASS_SW,rinfo.MASS_SED,rinfo.MASS_SED_DEEP,rinfo.PEC_SW,rinfo.PEC_SED]
                
                db_pars= ["MASS_SW","MASS_SED","MASS_SED_DEEP", "PEC_SW","PEC_SED"]
                db_pars = [par for par,hasPar in zip(db_pars,has_par) if hasPar]
                

                db_funcs=[np.take,np.take,np.take,np.mean,np.mean,np.mean,np.mean]
                db_funcs = [func for func,hasPar in zip(db_funcs,has_par) if hasPar]
                
                db_kwargs=[{'indices':-1,"axis":1},{'indices':-1,"axis":1},
                 {'indices':-1,"axis":1},{"axis":1},{"axis":1},
                 {"axis":1},{"axis":1}]
                db_kwargs = [kwargs for kwargs,hasPar in zip(db_kwargs,has_par) if hasPar]
     
                #TODO:
                db_pars = ["MASS_SW","MASS_SED","MASS_SED_DEEP","PEC_SW",
                           "PEC_SED","DEGR_SW","DEGR_SED"]
                db_funcs= [np.take,np.take,np.take, np.mean,np.mean,np.sum,np.sum]
                db_kwargs=[{'indices':-1,"axis":0},{'indices':-1,"axis":0}, 
                           {'indices':-1,"axis":0},{"axis":0},{"axis":0}, 
                           {"axis":0},{"axis":0}]
                                 
                sRN = STEPSRiverNetwork(rinfo.fpath,rinfo.key,rinfo.database,
                                        tstart=rinfo.begin,tend=rinfo.end,
                                        t_input_start=rinfo.t_input_start,
                                        timestep_input=rinfo.timestep_input,
                                        threshold_sw=rinfo.threshold_sw,
                                        threshold_sed=rinfo.threshold_sed,
                                        db_pars = db_pars,
                                        db_funcs = db_funcs,
                                        db_kwargs=db_kwargs,
                                         DEGHL_SW_0=cmp.DT50sw,
                                         DEGHL_SED_0=cmp.DT50sed,
                                         KOC=cmp.KOC,
                                         Temp0=cmp.Temp0,
                                         Q10=cmp.Q10,
                                         DENS=DENS,
                                         POROSITY=POROSITY,
                                         OC=OC,
                                         DEPTH_SED=DEPTH_SED,
                                         DEPTH_SED_DEEP=DEPTH_SED_DEEP,
                                         DIFF_L_SED=0.005,
                                         DIFF_L=0.005,
                                         DIFF_W=4.3*(10**-5),
                                         convertDeltaT=24.*60) # convert daily degrdataion value to minutes
                self.runs.append(sRN)
                
    def run(self,printres=False):  
        """
        Make a single model run with 'name' or comptue all runs listed in
        run List.
        """
        
        # make modle runs
        for run,rinfo in zip(self.runs,self.runlist):
            if rinfo.simulation == True:
                # make simulation
                run(printres)
                
                if rinfo.write_summary:
                    # write to disk
                    print("save reach-file")
                    run.write_reach_file(withHydro=rinfo.withHydro,
                                         agg=rinfo.aggregation_output,
                                         rule=rinfo.timestep_output)
            
                # close db
#                run.close_db()

    def preprocessing(self):
        """
        """
        pass

    def postprocessing(self):
        """
        """
        pass
    
    def message(self,s):
        """
        Writes a user message.
        """
        print(s)