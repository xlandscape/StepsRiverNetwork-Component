# -*- coding: utf-8 -*-
"""
.. module:: main
   :synopsis: Main script for running StepsRiverNetwork.
.. moduleauthor:: Sebastian Multsch <smultschw@knoell.com>

"""

import os
import argparse
from pathlib import Path
from RunFactory import RunFactory
import numpy as np
import h5py
import pandas as pd
from datetime import datetime
from Steps1234 import Steps1234
 
import timeit



def set_flows(fpath):
    dat = pd.read_csv(os.path.join(fpath,"HydroList_short.csv"))
    dat.set_index("key",inplace=True)
    dat.volume=10
    dat.flow=10
    dat.area=10
#    dat = dat.iloc[:24*365]
    dat.to_csv(os.path.join(fpath,"HydroList.csv"))

def runningMeanFast(x, N):
    return np.convolve(x, np.ones((N,))/N)[(N-1):]

def get_rolling_window(a, window):
    """
    Get a rolling window for a certain period.
    """
    shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
    strides = a.strides + (a.strides[-1],)
    return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)

def get_TWA(a,window,func_TWA):
    """
    Calculates time weighted average for given period.
    """
    
    # get TWA value for all reaches and return max
    reaches_twa = [np.mean(get_rolling_window(reach,window),axis=1).max() for reach in a]
        
    return func_TWA(reaches_twa)

def process_output(run,runtime,TWA=7,func_TWA=np.nanmax):
    """
    """
    

    
    fpath = run.fpath
    fname = run.key
    
    fh5 = h5py.File(os.path.join(fpath,fname+"_reaches.h5"), "r") 
    reachnames = run.reachnames

    
    s=""
    s+="#################################################################"+"\n"
    s+="# GENERAL"+"\n"
    s+="# PATH: " + fpath+"\n"
    s+="# NAME: " + fname +"\n"
    s+="# DATE: " + str(datetime.now())+"\n"
    s+="# RUNTIME: " + str(runtime)+"\n"
    s += "#\n"

    s+="#################################################################"+"\n"
    s+="# SCENARIO"+"\n"
    s+="# time start: " + str(run.t0)+"\n"
    s+="# time end: " + str(run.tn)+"\n"
    s+="# timestep hydrology: " + str(run.timestep_input)+"\n"
    s+="# timestep output: " + str(run.dt)+"\n"
    s+="# Number of reaches: " + str(len(reachnames))+"\n"
    s+="# River length: %.2f"%(sum(run.length))+"\n"
    s+="# Reach length min: %.2f"%(min(run.length))+"\n"
    s+="# River length max: %.2f"%(max(run.length))+"\n"
    s+="# River length mean: %.2f"%(np.mean(run.length))+"\n"
    s += "#\n"
    
    s+="#################################################################"+"\n"
    s+="# CATCHMENT TOTAL"+ "\n"
    
    t0 = run.t0
    tn = datetime(run.tn.year,run.tn.month,run.tn.day,run.tn.hour)
    times = list(pd.date_range(run.t0,run.tn,freq="1H"))
    i_t0 = times.index(t0)
    i_tn = times.index(tn) 
    i_t0_inp = run.time_input.index(t0)
    i_tn_inp = run.time_input.index(tn)    
    INPUT_SW = np.nansum(fh5["input_sw"][i_t0_inp:i_tn_inp,:])
    MASS_SW = np.nansum(fh5["MASS_SW"][i_tn,:])
    MASS_SED = np.nansum(fh5["MASS_SED"][i_tn,:]) 
    PEC_SW = np.nanmax(fh5["PEC_SW"][i_t0:i_tn,:]) 
    PEC_SED = np.nanmax(fh5["PEC_SED"][i_t0:i_tn,:])
    PEC_SW_TWA =  get_TWA(fh5["PEC_SW"][i_t0:i_tn,:].T,TWA,func_TWA)
    PEC_SED_TWA = get_TWA(fh5["PEC_SED"][i_t0:i_tn,:].T,TWA,func_TWA)  
    DEGR_SW = np.nansum(fh5["DEGR_SW"][i_t0:i_tn,:]) 
    DEGR_SED = np.nansum(fh5["DEGR_SED"][i_t0:i_tn,:]) 
    OUTFLOW = INPUT_SW - MASS_SW-MASS_SED-DEGR_SW-DEGR_SED
    s += "# %8s%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s"%("year","INPUT_SW","MASS_SW",
                                             "MASS_SED",
                                             "DEGR_SW","DEGR_SED","OUTFLOW","PEC_SW",
                                             "PEC_SED","PEC_SW","PEC_SED") + "\n"
    s += "# %8s%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s"%(" ","sum","end","end","sum","sum",
                                             "sum","max","max","7dTWA","7dTWA") + "\n"
    s += "# %8s%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s"%(" ","mg","mg","mg",
                                             "mg","mg","mg","ug/L","ug/L","ug/L","ug/L") + "\n"
    s += "%10s%10.2f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f"%("all",INPUT_SW,MASS_SW,
                                                           MASS_SED,DEGR_SW,
                                                           DEGR_SED,OUTFLOW,PEC_SW,
                                                           PEC_SED,PEC_SW_TWA,
                                                           PEC_SED_TWA)
    s += "\n"
    s += "#\n"
    
    s += "# %8s%10s%10s%10s%10s%10s%10s"%("year","INPUT_SW","MASS_SW",
                                             "MASS_SED",
                                             "DEGR_SW","DEGR_SED","OUTFLOW") + "\n"
    s += "# %8s%10s%10s%10s%10s%10s%10s"%(" ","sum","end","end","sum","sum",
                                             "sum") + "\n"
    s += "# %8s%10s%10s%10s%10s%10s%10s"%(" ","%","%","%",
                                             "%","%","%") + "\n"
    s += "%10s%10.2f%10.4f%10.4f%10.4f%10.4f%10.4f"%("all",100,MASS_SW/INPUT_SW*100.,
                                                           MASS_SED/INPUT_SW*100.,DEGR_SW/INPUT_SW*100.,
                                                           DEGR_SED/INPUT_SW*100.,OUTFLOW/INPUT_SW*100.)    
    s += "\n"
    s += "#\n"
        
    print(s)
        
    s+="#################################################################"+"\n"
    s+="# CATCHMENT ANNUAL"+ "\n"

    # write table header
    s += "# %8s%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s"%("year","INPUT_SW","MASS_SW",
                                             "MASS_SED",
                                             "DEGR_SW","DEGR_SED","OUTFLOW","PEC_SW",
                                             "PEC_SED","PEC_SW","PEC_SED") + "\n"
    s += "# %8s%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s"%(" ","sum","end","end","sum","sum",
                                             "sum","max","max","7dTWA","7dTWA") + "\n"
    s += "# %8s%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s"%(" ","mg","mg","mg",
                                             "mg","mg","mg","ug/L","ug/L","ug/L","ug/L") + "\n"

    # get reproting years
    times = list(pd.date_range(run.t0,run.tn,freq="1H"))
    years = np.unique([t.year for t in times])
    time_indices = []
    for y in years:
        time = [t for t in times if t.year == y]       
        
        if time[0] != time[-1]:
        
            time_indices.append([times.index(time[0]),times.index(time[-1])])
    time_indices_input = []
    for y in years:
        
        if y == run.t0.year:
            time = [t for t in run.time_input if t.year == y and t >= run.t0]  
        
        elif y == run.tn.year:
            time = [t for t in run.time_input if t.year == y and t <= run.tn]    
        else:        
            time = [t for t in run.time_input if t.year == y]  
        
        
        time_indices_input.append([run.time_input.index(time[0]),
                                   run.time_input.index(time[-1])])

    # create annual reports
    time_indices = [i for i in map(list, zip(*time_indices))] # transpose list
    time_indices_input = [i for i in map(list, zip(*time_indices_input))] # transpose list
    for year,istart,iend,inp_istart,inp_iend in zip(years,time_indices[0],
                                                    time_indices[1],
                                                    time_indices_input[0],
                                                    time_indices_input[1]):

        # get data
        INPUT_SW = np.nansum(fh5["input_sw"][inp_istart:inp_iend,:])
        MASS_SW = np.nansum(fh5["MASS_SW"][iend,:])
        MASS_SED = np.nansum(fh5["MASS_SED"][iend,:]) 
        PEC_SW = np.nanmax(fh5["PEC_SW"][istart:iend,:]) 
        PEC_SED = np.nanmax(fh5["PEC_SED"][istart:iend,:])
        PEC_SW_TWA =  get_TWA(fh5["PEC_SW"][istart:iend,:].T,TWA,func_TWA)
        PEC_SED_TWA = get_TWA(fh5["PEC_SED"][istart:iend,:].T,TWA,func_TWA)
        DEGR_SW = np.nansum(fh5["DEGR_SW"][istart:iend,:]) 
        DEGR_SED = np.nansum(fh5["DEGR_SED"][istart:iend,:]) 
        OUTFLOW = INPUT_SW-MASS_SW-MASS_SED-DEGR_SW-DEGR_SED

        # write annual data
        s += "%10i%10.2f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f"%(year,INPUT_SW,MASS_SW,
                                                           MASS_SED,DEGR_SW,
                                                           DEGR_SED,OUTFLOW,PEC_SW,
                                                           PEC_SED,PEC_SW_TWA,
                                                           PEC_SED_TWA)
        s += "\n"
    s += "#\n"

    s+="#################################################################"+"\n"
    s+="# REACHES"+ "\n"
    t0 = run.t0
    tn = datetime(run.tn.year,run.tn.month,run.tn.day,run.tn.hour)
    times = list(pd.date_range(run.t0,run.tn,freq="1H"))
    i_t0 = times.index(t0)
    i_tn = times.index(tn) 
    i_t0_inp = run.time_input.index(t0)
    i_tn_inp = run.time_input.index(tn)    
    
    INPUT_SW = np.nansum(fh5["input_sw"][i_t0_inp:i_tn_inp,:],axis=0)
    MASS_SW = fh5["MASS_SW"][i_tn,:]
    MASS_SED = fh5["MASS_SED"][i_tn,:] 
    PEC_SW = np.nanmax(fh5["PEC_SW"][i_t0:i_tn,:],axis=0) 
    PEC_SED = np.nanmax(fh5["PEC_SED"][i_t0:i_tn,:],axis=0)
    DEGR_SW = np.nansum(fh5["DEGR_SW"][i_t0:i_tn,:],axis=0) 
    DEGR_SED = np.nansum(fh5["DEGR_SED"][i_t0:i_tn,:],axis=0)   
    
    s += "# %8s%10s%10s%10s%10s%10s%10s%10s"%("key","INPUT_SW","MASS_SW",
                                             "MASS_SED",
                                             "DEGR_SW","DEGR_SED","PEC_SW","PEC_SED") + "\n"
    s += "# %8s%10s%10s%10s%10s%10s%10s%10s"%(" ","sum","end","end","sum","sum",
                                             "max","max") + "\n"
    s += "# %8s%10s%10s%10s%10s%10s%10s%10s"%(" ","mg","mg","mg",
                                             "mg","mg","ug/L","ug/L") + "\n"
    s+= "\n".join(["%10s%10.2f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f"%(r,inp,msw,msd,dsw,dsd,psw,psd) \
                   for r,inp,msw,msd,dsw,dsd,psw,psd in \
                   zip(reachnames,INPUT_SW,MASS_SW,MASS_SED,DEGR_SW,DEGR_SED,PEC_SW,PEC_SED)])
 
    # save file
    f = open(os.path.join(fpath,fname+".sum"),"w")
    f.write(s)
    f.close()        

    

    
    
    
def get_fdir(subfolder):
    """
    Return current working directoy and add subfolder to path.
    """
    fdir = os.path.abspath(os.path.join(
                    os.path.dirname(Path(__file__).parent),
                    *subfolder))
    return fdir

FLAGS = None

# todo:
# set value in sediment to last value
if __name__ == "__main__":
    
    start = datetime.now()
    
    # get command line arguments or use default
    parser = argparse.ArgumentParser()
    parser.add_argument('--folder',type=str, default=get_fdir(["projects"]),
                        help='Path of project folder.')
    parser.add_argument('--runlist',type=str,default='',
                        help='Name of model run')
    FLAGS, unparsed = parser.parse_known_args()


    # create run factory
    runfactory = RunFactory(FLAGS.folder,FLAGS.runlist)
     
    # setup model runs
    runfactory.setup()
    
    # conduct simulations
    runfactory.run(printres=False)
    
    # print overall processing time
    runtime = datetime.now()-start
    print("processing time",runtime)

    # process output
    print("process output")
    run = runfactory.runs[0]
    
    process_output(run=run,runtime=runtime)