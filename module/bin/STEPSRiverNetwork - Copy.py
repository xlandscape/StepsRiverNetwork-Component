# -*- coding: utf-8 -*-
"""
.. module:: STEPSRiverNetwork
   :synopsis: A river network version of Steps1234 (Klein, 2007).
.. moduleauthor:: Sebastian Multsch <smultschw@knoell.com>

"""

import h5py
import numpy as np
import pandas as pd
import os
from datetime import datetime
from Steps1234 import Steps1234

class STEPSRiverNetwork:
    """
    A river network version of Steps1234 (Klein, 2007).
    """
    def __init__(self,fpath:str,
                 key:str,
                 filetype:str="csv",
                 tstart:datetime = None,
                 tend:datetime = None,
                 t_input_start:datetime = None,
                 timestep_input:datetime = None,
                 threshold_sw:float=None,
                 threshold_sed:float=None,
                 db_pars:list= ["MASS_SW","MASS_SED","MASS_SED_DEEP","PEC_SW","PEC_SED"],
                 **kwargs):
        """
        :param fpath: Project folder
        :type fpath: str
        :param key: Project name
        :type key: str
        :param filetype: Filetype of input data (hdf or csv)
        :type filetype: str
        :param db_pars: List of parameters to written into database: 
            "MASS_SW","MASS_SED","MASS_SED_DEEP","PEC_SW","PEC_SED"
        :type db_pars: list[str]
        :returns: -
        :rtype: -
        """

        # set path variables
        self.fpath = os.path.join(fpath,key)
        self.key = key
        self.filetype = filetype
     
        self.t_input_start = t_input_start
        self.timestep_input = timestep_input
        
                
        
        self.threshold_sw = threshold_sw
        self.threshold_sed = threshold_sed
    

        # some user messages
        print("STEPSRiverNetwork simulation")
        print(self.fpath)
        print(self.key,"\n")

        # load reach and catchment list
        print("read input data")
        self.reachlist = pd.read_csv(os.path.join(self.fpath,"ReachList.csv"))
        self.reachnames = list(self.reachlist.key)
        self.n_reaches = len(self.reachnames)
        self.catchmentlist = pd.read_csv(os.path.join(self.fpath,"CatchmentList.csv"))

        
        # load hydrological input data
        fpath = os.path.join(self.fpath,"HydroList."+self.filetype)
        if self.filetype == "h5":
           
            # read hdf
            fh5 = h5py.File(fpath, 'r')
            self.flow = fh5["flow"][:]
            self.vol = fh5["volume"][:]
            self.area = fh5["area"][:]
            fh5.close()
            self.n_time_input = self.flow.shape[0]
    
        elif self.filetype == "csv":
            #read csv
            self.simulated = self.read_csv_to_pandas(fpath,
                                                     keys=["time","key"],
                                                     time_key="time",
                                                     time_format="%Y-%m-%dT%H:%M")
            self.n_time_input = len(self.simulated.index.levels[0])

            # reshape hydrological input data
            self.flow = self.simulated.flow.values.reshape([self.n_time_input,self.n_reaches])
            self.vol = self.simulated.volume.values.reshape([self.n_time_input,self.n_reaches])
            self.area = self.simulated.area.values.reshape([self.n_time_input,self.n_reaches])

        # set input dates
        self.time_input = list(pd.date_range(start=self.t_input_start,
                               periods=self.n_time_input,
                               freq=self.timestep_input))            
        
        # scale flow from m3/day to m3/min minutes  
        self.flow = np.divide(self.flow,1440.) 
        #TODO: use observed temperature          
        self.temp = np.ones([self.n_time_input,self.n_reaches],
                                dtype=np.float32) * 12 
        
        # get start time
        if tstart == None:
            self.t0_index =  0
            self.t0 =         self.time_input[0]
        else:
            self.t0 = tstart
            self.t0_index =  self.time_input.index(tstart)
        
        # get end time
        if tend == None:
            self.tn_index =  len(self.time_input)
            self.tn =         self.time_input[-1]  
        else:
            self.tn_index =  self.time_input.index(tend)+1
            self.tn =        tend
                
        self.dt =  self.time_input[1]-self.time_input[0]
        # if hourly, extent from 23:00 to 23:59
        if self.dt.seconds == 3600:
            self.tn = pd.Timestamp(year=self.tn.year,month=self.tn.month,
                                   day=self.tn.day,hour=self.tn.hour,
                                   minute=59)

        # crerate time index of output
        self.time = list(pd.date_range(start=self.t0,end=self.tn,freq="1Min"))
        self.n_time = len(self.time)
        

    
        # variable input solutes sw
        self.input_sw = np.zeros([self.n_time_input,self.n_reaches],dtype=np.float32)
    
        # spray drift
        self.spraydrift = pd.read_csv(os.path.join(self.fpath,"SprayDriftList.csv"))
        self.create_spraydrift()    
    
        # get geometry    
        rl_sort =  self.reachlist.sort_values("key")
        self.length = self.get_flowwidths()
        self.width = np.ones([self.n_reaches],dtype=np.float32) *rl_sort.bottomwidth.values

        # create efate model
        self.efate = Steps1234(self.length,self.width,**kwargs)    
    
        # set connections
        self.con = np.zeros([self.n_reaches,self.n_reaches],dtype=np.float32)
        self.create_cons()

        # create state variables   
        self.load_inflow = np.zeros([self.n_reaches],dtype=np.float32) 
        self.load_outflow = np.zeros([self.n_reaches],dtype=np.float32) 
        self.tb_mass = np.zeros([self.n_reaches],dtype=np.float32)     

        # create database
        self.db_pars = db_pars
        #self.db_dirs= [os.path.join(self.fpath,"steps_"+f+".npy") for f in self.db_pars]
        #self.db_dirs= [os.path.join(self.fpath,"steps_"+f+".h5") for f in self.db_pars]
        self.db=[]
        
#        # create db, delete existing files
#        for f in self.db_dirs: 
#            if os.path.isfile(f): 
#                os.remove(f)
#            npy = np.memmap(f,mode="w+",dtype=np.float32,
#                            shape=(self.n_time,self.n_reaches))
#            npy.flush()
#            self.db.append(npy)

        
        # create temproty files per hour
        dt_sub = int(self.dt.seconds / 60)
        for par in self.db_pars: 
            npy = np.zeros(dtype=np.float32,shape=(dt_sub,self.n_reaches))
            self.db.append(npy)


        # get datapath of db
        self.fname_h5 =  os.path.join(self.fpath,self.key+"_reaches.h5")
        
        # check if file exists and delete 
        if os.path.exists(self.fname_h5): 
            fh5 = h5py.File(self.fname_h5, 'r+')
            fh5.close()
            os.remove(self.fname_h5)
            
        # create hdf5 for final storage with a dset for each par 
        fh5 = h5py.File(self.fname_h5, 'w')
        for par in self.db_pars: 
            fh5.create_dataset(par, (self.tn_index-self.t0_index,self.n_reaches), 
                                 dtype="f4",compression="gzip",
                                 compression_opts=4)
        fh5.close()        

        # open db in append mode
        self.fh5 = h5py.File(self.fname_h5, "a") 
            
                    
    def printout(self,t:int,s:str,
                 vals:np.ndarray,
                 printres:bool=False):
        """
        Function to pritn user messages.
        :param t: Index of time step
        :type t: int
        :param s: Message
        :type s: str
        :param vals: Output values
        :type vals: np.ndarray
        :param printres: Printout detailed information
        :type printres: bool
        :returns: User message
        :rtype: str
        
        """
        if printres:
            return print(str(t)+ s + " ".join(["%.4f"%(i) for i in vals]))
        
    def __call__(self,printres:bool=False):
        """
        Make an efate simulation for a given set of hydrological data and 
        spray drift applications.
        
        mg/m3 = ug/L
        10**(-3)
        
        
        :param printres: Printout detailed information
        :type printres: bool
        :returns: -
        :rtype: -
        """



        print("simulate")
        start = datetime.now()
        tinternal = 0
        toutput = 0
        dt_sub = int(self.dt.seconds / 60)
        input_sw_at_zero = np.zeros(self.input_sw[0].shape,dtype=np.float32) 
        zero_state = False
        # run simulation
        for tinput in range(self.t0_index,self.tn_index,1):
            
            # get input
            input_sw_at_t = self.input_sw[tinput]
            
            # check condtions to skip sub-steps
            cond_inp = len(input_sw_at_t[input_sw_at_t>0])>0
            cond_sw = len(self.efate.PEC_SW[self.efate.PEC_SW>self.threshold_sw])>0
            cond_sed = len(self.efate.PEC_SED[self.efate.PEC_SED>self.threshold_sed])>0
            
            # run sub-ste
            for step in range(dt_sub):
                
                # run sub-set
                if (cond_inp or cond_sw or cond_sed):
                    
                    zero_state = False
        
                    if (tinternal % 1440 == 0):
                        print(self.time[tinternal],cond_inp,cond_sw,cond_sed,self.efate.PEC_SW.max(),self.efate.PEC_SED.max())
    
                    self.printout(tinput," input ",self.input_sw[tinput],printres)
                    
                    # 1. run steps1234
                    self.efate(input_sw_at_t,self.tb_mass,
                               self.temp[tinput],self.vol[tinput])
                    
                    #set input to zero
                    input_sw_at_t = input_sw_at_zero
                    
                    # 2. set mass 
                    self.tb_mass = self.efate.MASS_SW
                    self.printout(tinput," cmp start ",self.efate.MASS_SW,printres)
            
                    # 3. save results of current time step to memory
                    #self.save(tinternal,step)
                    self.save(toutput,step)                
                    
                    # 4. Calculate laod of reach outflow
                    self.load_outflow = self.efate.PEC_SW * self.flow[tinput]
                    self.printout(tinput," cmp out ",self.load_outflow,printres)
            
                    # 5. Substract outflow
                    self.tb_mass -= self.load_outflow
            
                    # 6. Sum up outflow loads as inflow for each reach according to connections 
                    self.load_inflow = np.sum(self.con * self.load_outflow,axis=1)
                    
                    # 7. Add inflow
                    self.tb_mass += self.load_inflow
                    self.printout(tinput," cmp next ",self.tb_mass,printres)
                    
                else:
                    if zero_state == False:
                        zero_state = True
                        for var in ["MASS_SW","MASS_SED","MASS_SED_DEEP","PEC_SW","PEC_SED"]:
                            vars(self.efate)[var][:] = 0. 
                        
                        
                        
                        
                
                
                
                # increase internal timestep
                tinternal += 1
                
            # increase output timestep
            toutput +=1

        # print runtime
        print("runtime",datetime.now()-start)  
        # final flush
        self.fh5.close()
 

        
    def get_flowwidths(self):
        """
        Calculates flow widths
        
        :returns: Array with flowwidth of each river segment.
        :rtype: np.array
        """     
        rl_sort =  self.reachlist.sort_values("key")
        outlet=self.catchmentlist[self.catchmentlist.component=="Outlet"]

        flowwidth = []
        for i in range(len(rl_sort)):         
            reach = rl_sort.iloc[i]
            
            # get coords
            x1 = reach.x
            y1 = reach.y
            if reach.downstream == "Outlet":
                x2 = outlet.x.values[0]
                y2= outlet.y.values[0]
            else:
                x2 = rl_sort[rl_sort["key"]==reach.downstream].x.values[0]
                y2 = rl_sort[rl_sort["key"]==reach.downstream].y.values[0]
                #calculate flow width               
            flowwidth.append(np.sqrt((x2-x1)**2 + (y2-y1)**2))
        return np.array(flowwidth,dtype=np.float32)
        
    def create_spraydrift(self):
        """
        Creates spray drift table by using the information from the SprayDriftList
        and the wetarea from hydrological data.
        
        :returns: -
        :rtype: -
        """   
        
        # set connections
        for i in range(len(self.spraydrift)):    
            event = self.spraydrift.iloc[i]
            # get reach index
            rind = self.reachnames.index(event.key)
            # get time index
            tind = self.time_input.index(pd.Timestamp(event.time))
            self.input_sw[tind,rind] = event.rate
        # multiply with wet area
        self.input_sw *= self.area    

    def create_cons(self):
        """
        Creates connetiosn between reaches according to ReachList.
        
        :returns: -
        :rtype: -
        """   
        # set connections
        for i,name in enumerate(self.reachnames):    
            # get inflow from others
            inflow = self.reachlist[self.reachlist.downstream==name].key.values
            # get index of inflows
            inflow_ind = [self.reachnames.index(reach) for reach in inflow]
            #set connection
            for ind in inflow_ind:
                self.con[i,ind]=1
        
#    def save(self,t):
#        """
#        Saves current timestep to memory and flushed data to disk each hourly
#        modelling timestep.
#        
#        :returns: -
#        :rtype: -
#        """   
#        for par,ds in zip(self.db_pars,self.db):
#            ds[t] = vars(self.efate)[par]
#            if t % 3600 == 0:
#                for ds in self.db:
#                    ds.flush()

    def save(self,texternal,step):
        """
        Saves current timestep to memory and flushed data to disk each hourly
        modelling timestep.
        
        :returns: -
        :rtype: -
        """   
       
        # save all required  parameter
        for par,ds in zip(self.db_pars,self.db):
            # save actual step per minute into tempory numpy-array
            ds[step] = vars(self.efate)[par]
            
        # calc hourly summary and reset tempory numpy-array
        if step == 59: 
#            f= h5py.File(self.fname_h5, "a")    
            for par,ds in zip(self.db_pars,self.db):
                dset = self.fh5[par]
                dset[texternal] = ds.max(axis=0)
                ds*=0 # set array to zero
#            f.close()
                
    def write_reach_file(self,withHydro=False,agg=None,rule="1H"):
        """
        Create pandas dataframe from database using the input hydrological data
        and the simulated efate output. The simulation data per minute is 
        aggregated to hours by calculating min,max or mean per hour.
        
        :returns: Datatframe wit hhydrology and efate data.
        :rtype: pandas.DataFrame
        """
        
        # create list with keys
        keys = np.array([(t,r) for r in self.reachnames for t in self.time])

        # create dataframe
        df = pd.DataFrame()
        df["key"] = keys[:,1]
        df["time"] = keys[:,0]
        for par,ds in zip(self.db_pars,self.db):
            df[par] = ds.T.flatten()

        # set indices
        df.set_index(["time","key"],inplace=True)
        
        # aggregate
        if agg!=None:
            gb=df.groupby(level=["key"])
            
            if agg == "min":
                df=gb.apply(lambda x: x.reset_index(level=["key"]).resample(rule=rule).pad())

            elif agg == "max":
                df=gb.apply(lambda x: x.reset_index(level=["key"]).resample(rule=rule).max())

            elif agg == "mean":
                df=gb.apply(lambda x: x.reset_index(level=["key"]).resample(rule=rule).mean())
            
            df=df.swaplevel()[["MASS_SW","MASS_SED","MASS_SED_DEEP","PEC_SW","PEC_SED"]]
            agg = "_" + agg

        
        # sort data
        df.sort_index(level=[0], ascending=True,inplace=True)
         
        # join with hydro
        if withHydro:
            print("merge")
            df = self.simulated.merge(df,left_index=True,right_index=True)
        
        if agg==None:
            agg=""
            
        # save to disk
        df.to_csv(os.path.join(self.fpath,"steps_reaches"+agg+".csv"),
                  date_format="%Y-%m-%dT%H:%M")
        
        return df

    def read_hdf_to_pandas(self,fname="",keys=None, 
                           dset=None,
                           time_key=None,
                           time_format="%Y-%m-%dT%H:%M",
                           convert_byte=None):
        """
        Reads a hdf5 file with h5py and returns pandas dataframe.

        By defining a keys and/or a time_key and index is established and
        string-dates are converted to datetime.
        convert_byte is a variable for a workaround of a datatype issue between
        Python3 numpy and h5py.

        :param fname: Path and name of file with extention.
        :type fname: str
        :param keys: List with keys.
        :type keys: list(str)        
        :param dset: Name of dataset in hdf5 file.
        :type keys: str  
        :param time_key: Name of time key.
        :type time_key: str        
        :param time_format: Datetime format, e.g."%Y-%m-%dT%H:%M".
        :type time_format: str        
        :param convert_byte: Variables where fomrat must be converted..
        :type convert_byte: list(str)   
        
        :returns: Data frame
        :rtype: pandas.DataFrame        
        """
        # read hdf and extract all data into numpy array
        npy = h5py.File(fname, 'r')[dset][:]
        # create DataFrame from numpy arra
        data = pd.DataFrame(npy)
        # convert byte-string to strings @TODO: slow
        for par in convert_byte:
            data[par]= data[par].apply(lambda x: str(x,'utf-8'))
        # convert time string to datetime objects
        if time_key != None:
            data ["time"]=pd.to_datetime(data [time_key],format=time_format)
        # set keys
        if keys != None:
            data.set_index(keys,inplace=True)
        return data
        
    def read_csv_to_pandas(self,fname="",keys=None,time_key=None,
                           time_format="%Y-%m-%dT%H:%M",**kwargs):
        
        """
        Reads a csv file with h5py and returns pandas dataframe.
    
        By defining a keys and/or a time_key and index is established and
        string-dates are converted to datetime.
        convert_byte is a variable for a workaround of a datatype issue between
        Python3 numpy and h5py.
    
        :param fname: Path and name of file with extention.
        :type fname: str
        :param keys: List with keys.
        :type keys: list(str)        
        :param time_key: Name of time key.
        :type time_key: str        
        :param time_format: Datetime format, e.g."%Y-%m-%dT%H:%M".
        :type time_format: str        
        
        :returns: Data frame
        :rtype: pandas.DataFrame        
        """
        # read file
        data = pd.read_csv(fname,**kwargs)
        # create dattime from string
        if time_key != None:
            data ["time"]=pd.to_datetime(data [time_key],format=time_format)
        # set keys
        if keys != None:
            data.set_index(keys,inplace=True)
        return data    