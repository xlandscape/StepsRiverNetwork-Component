# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 11:33:24 2019

@author: smh
"""

# python native
import os
from calendar import monthrange
from datetime import timedelta
import shutil
import copy
import math

# other
import pandas as pd
import numpy as np

# plotting
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.axisartist as AA
from mpl_toolkits.axes_grid1 import host_subplot
from mpl_toolkits.axes_grid.parasite_axes import SubplotHost
from matplotlib import collections  as mc
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
        

# own
from .Parameter import ParameterList
from .Tools import TimeSeries,IO,Plots
from .Signatures import Signatures

# define Python user-defined exceptions
class Error(Exception):
   """Base class for other exceptions"""
   pass

class MissingFile(Error):
   """Raised when a file is missing"""
   pass

class PostProcessing(TimeSeries,IO,Plots,Signatures):
    
    def __init__(self,fpath,
                 fname,
                 time_format="%Y-%m-%dT%H:%M"):
        
        
        # implement parent classes
        TimeSeries.__init__(self)
        IO.__init__(self)
        Plots.__init__(self)
        Signatures.__init__(self)
        
        # set catchment project
        print("postprocessing:", os.path.join(fpath) )
        self.fpath = fpath
        self.fname = fname

        self.time_format = time_format

        # load observed data
        fpath_obs =  os.path.join(self.fpath,"observation.csv")
        if os.path.isfile(fpath_obs):
            #print("load", fpath_obs)
            self.observed = pd.read_csv(fpath_obs)
            self.observed ["time"]=pd.to_datetime(self.observed ["time"],
                          format=self.time_format)
            self.observed.set_index(["time","key"],inplace=True)
                    
            # get gauging station
            self.gaugingstation = self.observed.index.levels[1][0]
        else:
            #print("file not exists",fpath_obs)
            self.observed = None
            self.gaugingstation = None
        
        # load spraydrift list
        self.spraydrift = pd.read_csv(os.path.join(self.fpath,"SprayDriftList.csv"))
        
        # load reaches
        self.reaches = pd.read_csv(os.path.join(self.fpath,"ReachList.csv"))
        self.reaches.set_index("key",inplace=True)
        
        # laod catchmentlist
        self.catchmentlist = pd.read_csv(os.path.join(self.fpath,"CatchmentList.csv"))
        self.catchmentlist.set_index("key",inplace=True)        
    
        try:
            # get filetype
            if os.path.isfile(os.path.join(self.fpath,self.fname+"_reaches.csv")):
                self.filetype = "csv"
                fpath_sim = os.path.join(self.fpath,self.fname + "_reaches."
                                         +self.filetype)
                self.simulated = self.read_csv_to_pandas(fname=fpath_sim,
                                             keys=["time","key"],
                                             time_key="time",
                                             time_format="%Y-%m-%dT%H:%M")
            elif os.path.isfile(os.path.join(self.fpath,self.fname+"_reaches.hdf")):
                self.filetype = "hdf"
                fpath_sim = os.path.join(self.fpath,self.fname + "_reaches."
                             +self.filetype)
                self.simulated =self.read_hdf_to_pandas(fname=fpath_sim,
                                            keys=["time","key"], 
                                            dset="reaches",
                                            time_key="time",
                                            time_format="%Y-%m-%dT%H:%M",
                                            convert_byte=["time","key"])
            else:
                raise MissingFile
        except MissingFile:
            print("MissingFile: " + os.path.join(self.fpath,self.fname+"_reaches"))
            

    


       
            
#        # climfile
#        self.clim_file = pd.read_csv(os.path.join(self.fpath,"%s.csv" % ([row.key for row in self.catchment.inpData.ClimateList][0])))
#        self.clim_file["time"]=pd.to_datetime(self.clim_file["key"],format=self.time_format)
#        self.clim_file.set_index(["time"],inplace=True) 
#
#
#    def calc_signatures(self):
#        
#        results_path = os.path.join(self.fpath,"Signatures")
#        Signatures.createFolder(self, results_path) 
#
#        try:
#            # checks if observation period covers simulation period (years)
#            pct = Signatures.check_data_inputs(self, self.simulated, self.observed)
#            
#            if pct == 100:            
#                results_path = os.path.join(self.fpath,"Signatures")
#                Signatures.createFolder(self, results_path) 
#                print('Folder Signatures created')
#                Signatures.calc_signatures(self, self.simulated, self.observed, self.gaugingstation, self.fname, results_path, self.clim_file)
#                print('Signatures for complete timeseries calculated')
#                Signatures.calc_yearly_TS(self, self.simulated, self.observed, self.gaugingstation, self.fname, results_path)
#                print('Signatures for yearly timeseries calculated')
#                print('### DONE ###')
#            else:
#                print('simulation period not in observation period')
#                Signatures.calc_signatures_sim_only(self, self.simulated, self.gaugingstation, self.fname, results_path, self.clim_file)
#        except:
#            Signatures.calc_signatures_sim_only(self, self.simulated, self.gaugingstation, self.fname, results_path, self.clim_file)
#            

    def get_zero_flow_reaches(self,stats=None):
        """
        Get all reaches with at least one value with flow == 0 and make a plot.
        """
        
        #######################################################################
        # process stats
        
         # get list with reaches with zero flow
        reaches_no_water = stats[stats.flow["min"]==0].index.tolist()
        
        # select data with zero flow
        vals = self.select_by_indices(params = ["depth","flow"],
                              times=slice(None),
                              keys=reaches_no_water)
    
        # group reaches with zero flow
        gb = vals["flow"].groupby(level=1)
        groups = dict(list(gb))
        times = vals.index.levels[0].to_pydatetime()

        # merge stats with coordiantes of reaches
        res = stats.merge(self.reaches, 
                            left_index=True, right_index=True,
                            suffixes=('', ''))
        
        # merge with downstream reach coordiantes
        res = res.merge(self.reaches, 
                            left_on='downstream', 
                            right_index=True,
                            suffixes=('_upstream', '_downstream'))    
        
        #######################################################################
        # plot time series
        
        # plot time series
        fig, ax = plt.subplots(figsize=(15,7))
        for group in groups:
            ax.plot(times,groups[group].values)
        ax.grid(True)
        ax.set_xlabel("time")    
        ax.set_ylabel("Flow [m$^3$ sec$^{-1}$]")
        fname = os.path.join(self.fpath,"reaches_no_water.png")
        fig.savefig(fname,dpi=300)
        fig.clf()
      
        #######################################################################
        # plot map
        
        x_upstream = res["x_upstream"].values
        y_upstream= res["y_upstream"].values
        vals=res[('flow', 'min')].values
        vals[vals>0]=1
        vmin=vals.min()
        vmax=vals.max()
           
        # make figure
        fig= plt.figure(figsize=(15,15))
        ax1 = fig.add_axes([0.1,0.1,0.8,0.8]) # x,y, lenght, height
    
        #create a colorbar
        cmap = matplotlib.cm.get_cmap('jet_r')
        norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
        
        # create map
        ax1.set_xlim((x_upstream.min()),x_upstream.max())
        ax1.set_ylim((y_upstream.min(),y_upstream.max()))
        ax1.set_xlabel("X [m]")
        ax1.set_ylabel("Y [m]")
        ax1.grid(True)    
    
        # add data
        lines = [[x,y] for x,y in zip(res[["x_upstream","y_upstream"]].values,
                                     res[["x_downstream","y_downstream"]].values)]
        colors = [cmap(norm(val)) for val in vals]
        lc = mc.LineCollection(lines, colors=colors, linewidths=2)
        ax1.add_collection(lc)
        ax1.autoscale()
        ax1.margins(0.1)
    
        # add key of reaches with zero values
        text=res[res[('flow', 'min')]==0]
        xvals = (text["x_upstream"].values+text["x_downstream"].values) / 2.
        yvals = (text["y_upstream"].values+text["y_downstream"].values) / 2.
        svals = text.index
        for x,y,s in zip(xvals,yvals,svals):
            ax1.text(x,y,s,verticalalignment='bottom', 
                     horizontalalignment='right',
                            color='r', fontsize=10, 
                            bbox=dict(facecolor='r', edgecolor='None', 
                                      boxstyle='round,pad=.2',alpha=.25))
        # make colormap
        divider = make_axes_locatable(ax1)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm._A = []
        sm
        cb=fig.colorbar(sm,cax=cax, orientation='vertical')
        cb.set_label("0: flow=0 flow; 1=flow>0")
    
        plt.close("all")
        fig.savefig( os.path.join(self.fpath,"map_no_water.png"),dpi=300,transparent=True) 
    
        fig.clf()    
    

    def plot3D_params(self):
        pass

    
        #
        #    # make more specifci post processing
        #    pstPrc = LandscapeModel.utils.PostProcessing(os.path.join(FLAGS.fdir,FLAGS.key),
        #                                                         FLAGS.key,time_format="%Y-%m-%dT%H:%M")
        
        #
        #    # calculate stats and save file
        #    stats = pstPrc.get_stats(stats= ['mean', 'median',"min","max"], 
        #                             params=['depth', 'volume', 'flow', 'area','artificialflux'])
        #    # merge stats with coordiantes of reaches
        #    res = stats.merge(pstPrc.reaches,  left_index=True, right_index=True, suffixes=('', ''))
        #    # merge with downstream reach coordiantes
        #    res = res.merge(pstPrc.reaches, 
        #                        left_on='downstream', 
        #                        right_index=True,
        #                        suffixes=('_upstream', '_downstream'))    
        #    # merge with slope and flowdsith
        #    flowwidths = pstPrc.get_flowwidths()
        #    slopes = pstPrc.get_slopes()
        #    res = res.merge(flowwidths, 
        #                        left_index=True, 
        #                        right_index=True,
        #                        suffixes=('_fw', '_fw'))       
        #    res = res.merge(slopes, 
        #                        left_index=True, 
        #                        right_index=True,
        #                        suffixes=('_sl', '_sl'))       
        #
        #    font = {'family' : 'normal',
        #        'weight' : 'normal',
        #        'size'   : 10}
        #    matplotlib.rc('font', **font) 
        #    
        #
        #    slope =   res.slope*-1
        #    volume = res["volume","max"].values
        #    flow = res["flow","max"].values/ 86400
        #    depth = res["depth","max"].values 
        #
        #    # plot time series
        #    fig = plt.figure(figsize=(7,5))
        #    ax = fig.add_subplot(111, projection='3d')
        #    ax.scatter(slope,volume,flow, s=50,edgecolor="None", marker="o",c=depth,
        #                          alpha=0.5,cmap="jet_r",lw=0)
        #    cbar=plt.colorbar(sc)
        #    ax.set_xlabel("Slope [-]")
        #    ax.set_ylabel("Volume [m³]")
        #    ax.set_zlabel("Flow [m³/sec]")
        #    cbar.set_label("Depth [m]")
        #    ax.grid(True)
        #
        ###
        ##    dat = res[(slope<0.00002)&(volume>600)]
        ##
        ##    xs= (dat.slope*-1).values
        ##    ys=dat["volume","max"].values
        ##    zs=dat["flow","max"].values/ 86400
        ##    ls=dat.index.values
        ##    for x,y,z,l in zip(xs,ys,zs,ls):
        ##        ax.text(x,y , z , l,"y")
        #
        #
        #    plt.close("all")
        #    fig.savefig( os.path.join(pstPrc.fpath,"00_slopes2.png"),
        #        dpi=300,transparent=True) 
        #   
    def map_stats(self,stats=None):
        """
        """

        # calculate stats and save file
        stats = self.get_stats(stats= ['mean', 'median',"min","max"],
                      params=['depth', 'volume', 'flow', 'area'])
        
        # merge stats with coordiantes of reaches
        res = stats.merge(self.reaches, 
                            left_index=True, right_index=True,
                            suffixes=('', ''))
        
        # merge with downstream reach coordiantes
        res = res.merge(self.reaches, 
                            left_on='downstream', 
                            right_index=True,
                            suffixes=('_upstream', '_downstream'))    
    

        # make map of all stats
        for param in stats.columns.levels[0]:
            for statval in stats.columns.levels[1]:
                
                print("plot",param,statval)
          
                ###############################################################
                # plot map
                
                # get values 
                vals=res[(param, statval)].values

                   
                # make figure
                fig= plt.figure(figsize=(15,15))
                
                
                ax1 = fig.add_axes([0.1,0.1,0.8,0.8]) # x,y, lenght, height
            
                #create a colorbar
                cmap = matplotlib.cm.get_cmap('jet_r')
                norm = matplotlib.colors.Normalize(vmin=vals.min(), vmax=vals.max())
                
                # set axis propterties
                xmin = res[["x_upstream","x_downstream"]].min().min()
                xmax=res[["x_upstream","x_downstream"]].max().max()
                ymin = res[["y_upstream","y_downstream"]].min().min()
                ymax=res[["y_upstream","y_downstream"]].max().max()    
                ax1.set_xlim(xmin,xmax)
                ax1.set_ylim(ymin,ymax)
                ax1.set_xlabel("X [m]")
                ax1.set_ylabel("Y [m]")
                ax1.grid(True)    
            
                # plot river segemnts
                lines = [[x,y] for x,y in zip(res[["x_upstream","y_upstream"]].values,
                                         res[["x_downstream","y_downstream"]].values)]
                colors = [cmap(norm(val)) for val in vals]
                lc = mc.LineCollection(lines, colors=colors, linewidths=2)
                ax1.add_collection(lc)
                ax1.autoscale()
                ax1.margins(0.1)
            
                # add key of reaches with zero values
                if statval == "min":
                
                    text=res[res[(param, statval)]<np.percentile(res[(param, statval)],1)]
                else:
                    text=res[res[(param, statval)]>np.percentile(res[(param, statval)],99)]
                    
                
                xvals = (text["x_upstream"].values+text["x_downstream"].values) / 2.
                yvals = (text["y_upstream"].values+text["y_downstream"].values) / 2.
                svals = text.index
                
                for x,y,s in zip(xvals,yvals,svals):
                    ax1.text(x,y,s,verticalalignment='bottom', 
                             horizontalalignment='right',
                                    color='k', fontsize=10, 
                                    bbox=dict(facecolor='k', edgecolor='None', 
                  
                                        boxstyle='round,pad=.2',alpha=.25))
                
            
                # make colormap
                divider = make_axes_locatable(ax1)
                cax = divider.append_axes('right', size='5%', pad=0.05)
                sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
                sm._A = []
                sm
                cb=fig.colorbar(sm,cax=cax, orientation='vertical')
                cb.set_label(" ".join([statval,param]),fontweight="bold",fontsize=12)
                
#                # north arrow
#                x_arrow = xmax*0.998
#                y_arrow = ymin * 1.008
#                ax1.arrow(x_arrow,y_arrow,0, 200, 
#                          head_width=800, 
#                          head_length=300, 
#                          shape ="right",
#                          overhang=5,
#                          fc='k', ec='k',
#                          alpha=1,label=ax1.text(x_arrow,y_arrow, "N", fontsize=20,
#                                                   fontweight="bold",alpha=1))
#            
                plt.close("all")
                fig.savefig( os.path.join(self.fpath,"map"+param+"_"+statval+".png"),
                    dpi=300,transparent=True) 
            
                fig.clf()    
                 

    def plot_percentile_over_time(self,params=[]):
    
        for param in params:
            time = self.simulated.index.get_level_values(0)
            gb = self.simulated[param].groupby([time.year,time.month,time.day])
            
            Q1 = gb.quantile(.25)
            Q2 = gb.quantile(.50)
            Q3 = gb.quantile(.75)
            IQR = Q3-Q1
            # The position of the whiskers is set by default to 1.5 * IQR (IQR = Q3 - Q1) 
            # from the edges of the box. Outlier points are those past the end of the whiskers.
            whiskers_high = Q3 + (IQR*1.5)
            whiskter_low =  Q1-(IQR*1.5)
    
            # plot time series
            fig, ax = plt.subplots(figsize=(10,7))
        
            # plot median, IQR and whiskers
            ax.plot(Q2.values,label="Median",color="b",linewidth=2)
            ax.plot(whiskers_high.values,label="Outlier",color="b",
                    linewidth=2,linestyle="--",alpha=.25)    
            plt.fill_between(np.arange(0,len(Q1),1),Q1,Q3, facecolor='b', 
                             edgecolors="None",alpha=.25)
        
        
            # plot statistics of fliers
            doy=0
            n = []
            for dat in gb:
                w = whiskers_high.iloc[doy]
                vals=dat[1][dat[1]>w]
                ax.plot(doy,vals.values.max(),marker="+",color="b",alpha=.25,linewidth=0)
                doy+=1
                n.append(len(vals.values))
                
            # format axis 
            ax.grid(True)
            ax.set_xlabel("Day of year")
            ax.set_ylabel(param)
    
            
            # legend 
            median = mlines.Line2D([],[],color='b', label="Median",alpha=1,linewidth=2,
                                   linestyle="-")
            iqr = mpatches.Patch(color='b', label='Interquartile range (IQR)'
                                 ,alpha=.25,linewidth=0)
            whiskers = mlines.Line2D([],[],color='b', label="Whiskers (>Q3+1.5*IQR)",
                                     alpha=.5,linewidth=2,linestyle="--")
            fliers = mlines.Line2D([],[],color='b', label="Fliers (max(>Whiskers))",
                                   alpha=.25,linewidth=0,marker="+")
#            lg=plt.legend(handles=[median,iqr,whiskers,fliers],ncol=1,loc="upper left",fontsize="large",
#                       title="Rummen",title_fontsize=12,facecolor="b",fancybox=True )
            
            
            lg=plt.legend(handles=[median,iqr,whiskers,fliers],ncol=2,fontsize="large",
                   title=None,title_fontsize=12,facecolor="b",fancybox=True,
                   bbox_to_anchor=(0.1, 1))
            
            lg.get_frame().set_alpha(0.1)  
            
                
            plt.close("all")
            fig.savefig( os.path.join(self.fpath,"ts_percentiles_"+param+".png"),
                dpi=300,transparent=True) 
        
            fig.clf()   


    def get_stats(self,stats= ['mean', 'median',"min","max"],
                  params=['depth', 'volume', 'flow', 'area']):
        """
        Calculates stats for a given set of paramter.
        """
        agg = dict([(par,stats) for par in params])
        stats = self.simulated.groupby(["key"]).agg(agg)
        stats.to_csv(os.path.join(self.fpath,"catchment_stats.csv"))
        return stats

    def select_by_indices(self,
                          params = [""],
                          times=slice(None),
                          keys=slice(None)):
        """
        Select a set of paramter by time and key index.
        """
        return self.simulated.loc[(times,keys),params]
    
    def performance(self,var):
        """
        makes a timeseries plot of observed and simulated values.
        """
        
        # get observed and simulated values
        obs = self.observed.xs(self.gaugingstation ,level=1)[var] 
        sim = self.simulated.xs(self.gaugingstation ,level=1)[var]  
        time = self.simulated.index.levels[0]
        
        # make plot
        fname = os.path.join(self.fpath,var+"_"+self.gaugingstation +".png")
        self.figure_performance(time,obs,sim,
                                "Flow [m$^3$ day$^{-1}]$","%b", fname)
    
    
    def catchment_hydrology(self):
        """
        Creates several histograms of hydrology related to compound application.
        """
        
        # get time
        time = self.simulated.index.get_level_values(0)
        key = self.simulated.index.get_level_values(1)
        
        # select values from reaches which receive drift at the day of application
        cond1 = ( key.isin(self.spraydrift.key.tolist()))
        cond2 = (time.day == 10) & (time.month == 5) & ( key.isin(self.spraydrift.key.tolist()))
        
        # get flowdiwidths of reaches
        flowwidth = self.get_flowwidths()
                
        # make join between flowwidth and simulations and normalise volume and area
        join = self.simulated.join(flowwidth,how="left")
        volume = join.volume / join.flowwith
        area = join.area / join.flowwith
        
        # make plot
        fig = plt.figure(figsize=(10, 8))
        
        # plot water depth
        ax = fig.add_subplot(221)
        ax.grid(True,alpha=.2,color="k")
        self.ax_hist(ax,data=[self.simulated.depth[cond1],self.simulated.depth[cond2]],
             colors=["k","orange"],labels=["Catchment","Appl. day"],
             xlabel="Depth [m]",text="(a)",alpha=.5,normed=True)
        
        # plot volume  
        ax = fig.add_subplot(222)
        ax.grid(True,alpha=.2,color="k")
        self.ax_hist(ax,data=[volume[cond1],volume[cond2]],
             colors=["k","orange"],labels=["Catchment","Appl. day"],
             xlabel="Volume [m$^3$ m$^{-1}$]",text="(b)",alpha=.5,normed=True)    
    
        # plot wet area  
        ax = fig.add_subplot(223)
        ax.grid(True,alpha=.2,color="k")
        self.ax_hist(ax,data=[area[cond1],area[cond2]],
             colors=["k","orange"],labels=["Catchment","Appl. day"],
             xlabel="Wet area [m$^2$ m$^{-1}$]",text="(c)",alpha=.5,normed=True)    
    
          # plot wet area  
        ax = fig.add_subplot(224)
        ax.grid(True,alpha=.2,color="k")
        self.ax_hist(ax,data=[self.simulated.flow[cond1],self.simulated.flow[cond2]],
             colors=["k","orange"],labels=["Catchment","Appl. day"],
             xlabel="Flow [m$^3$ day$^{-1}$]",text="(d)",alpha=.5,normed=True)  
        
        # save 
        fname = os.path.join(self.fpath,"catchment_hydrology.png")
        fig.savefig(fname,dpi=300)
        fig.clf()
    
    def catchment_efate(self,day,intervals,maxSW=1,maxSED=1):
        """
        Creates several histograms of efate related to compound application.
        """
        
        # get time
        time = self.simulated.index.get_level_values(0)
        key = self.simulated.index.get_level_values(1)
        
        conds = []
        for intv in intervals:
            conds += [(time==(day+timedelta(seconds=3600*intv))) & ( key.isin(self.spraydrift.key.tolist()))]
        
        # make plot
        fig = plt.figure(figsize=(10, 8))
        
        # create colors
        cmap = plt.cm.get_cmap("jet", len(intervals)+1)
        colors = [cmap(i) for i in range(len(intervals))]
        
        # create labels
        labels = ["app.time + %ih"%i for i in intervals]
        
        # plot water depth
        ax = fig.add_subplot(211)
        ax.grid(True,alpha=.2,color="k")

        self.ax_hist(ax,data= [self.simulated.PEC_SW[cond] for cond in conds],
             colors=colors,labels=labels,
             xlabel="PEC$_{SW}$ [$\mu$g L$^{-1}$]",text="(a)",alpha=.5,normed=True,
              histtype='step', cumulative=-1)
        ax.set_xlim(0, maxSW)
        
        # plot maximum vaulue
        ax.text(0.8, 0.98, "max:%.3f"%(self.simulated.PEC_SW.max()) + " $\mu$g L$^{-1}$", 
                  horizontalalignment='left', 
                  verticalalignment='top', transform=ax.transAxes,
                  fontsize=10,fontweight="bold")
        
        # plot volume  
        ax = fig.add_subplot(212)
        ax.grid(True,alpha=.2,color="k")
        self.ax_hist(ax,data= [self.simulated.PEC_SED[cond] for cond in conds],
             colors=colors,labels=labels,
             xlabel="PEC$_{SED}$ [$\mu$g L$^{-1}$]",text="(b)",alpha=.5,normed=True,
              histtype='step', cumulative=-1)  
        ax.set_xlim(0, maxSED)
        
        # plot maximum value
        ax.text(0.8, 0.98, "max:%.3f"%( self.simulated.PEC_SED.max()) + " $\mu$g L$^{-1}$", 
                  horizontalalignment='left', 
                  verticalalignment='top', transform=ax.transAxes,
                  fontsize=10,fontweight="bold")
   
        # save 
        fname = os.path.join(self.fpath,"catchment_efate.png")
        fig.savefig(fname,dpi=300)
        fig.clf()

    def reach_getTimeseries(self,reach,variable=None,tstart=None,tend=None):
        """
        Return the timeseries of a specific reach. tstart/tend are optional
        to set a time snippet.
        
        :param reach: Name of reach.
        :type reach: str
        :param variable: Name of variable.
        :type variable: str
        :param tstart: Start time.
        :type tstart: datetime
        :param tend: End time.
        :type tend: datetime  
        
        :returns: Timeseries
        :rtype: pandas.TimeSeries 
        """
        if tstart == None:
            dat =  self.simulated.xs(reach,level=1)
        else:
            dat = self.simulated.xs(reach,level=1).loc[tstart:tend]
        
        if variable != None:
            return dat[variable]
        else: return dat

    def branch_hydrology_and_efate(self,par1,lab1,par2,lab2,
                                   reach_start,reach_end,
                                   tstart,tend,tintverval,
                                   timescale="daily",
                                   vmin=None,vmax=None,fname=""):
        """
        Plot 3D bar chart with with reaches on x-axis, time on y-axis and par1
        on z-axis as well as par2 as color of the bars.
        
        :param par1: Parameter to be plotted on z-axis
        :type par1: float
        :param lab1: Axis label of first parameter.
        :type lab1: str
        :param par1: Parameter to be plotted as color of bars
        :type par1: float   
        :param par2: Axis label of second parameter.
        :type lab2: str
        :param reach_start: Reach at the start of branch.
        :type reach_start: str
        :param reach_end: Reach at the end of branch.
        :type reach_end: str    
        :param tstart: Start time.
        :type tstart: datetime
        :param tend: End time.
        :type tend: datetime  
        :param tintverval: Invervall of time axis.
        :type tintverval: int
        :param vmin: Minimum value colormap par2.
        :type vmin: float
        :param vmax:  Maximum value colormap par2.
        :type vmax: float            
        
        :returns: -
        :rtype: -
        """

        
        # find reaches between start and end
        reaches = [reach_start]                     
        complete = False
        while not complete:
            # get current reach
            reach = reaches[-1]
            # get downstream storage
            downstream = self.reaches.loc[reach].downstream # rList[reach][0].__dict__["downstream"]
            # stop at outlet
            if downstream == "Outlet":
                complete = True
            else:
                # stop at end reach
                if downstream == reach_end:
                    reaches.append(downstream)
                    complete = True
                # go on when not outlet and not end
                else:
                    # get next reach
                    reaches.append(downstream)

        # get flowwidth
        flowwidths = self.get_flowwidths()
        flowwidths = [flowwidths.loc[reach].flowwith for reach in reaches]
        flowwidths_cumulative = np.cumsum(flowwidths)
    
    	# get data snippet for timeperiod of par1 and par2
        dat1 = pd.DataFrame([self.simulated.xs(reach,level=1).loc[tstart:tend][par1].values for reach in reaches])
        dat2 = pd.DataFrame([self.simulated.xs(reach,level=1).loc[tstart:tend][par2].values for reach in reaches])


        # format dat1 for plotting
        xpos=np.arange(dat1.shape[0])
        ypos=np.arange(dat1.shape[1])
        yposM, xposM = np.meshgrid(ypos+0.5, xpos+0.5)
        zpos=np.zeros(dat1.shape).flatten()
        dx = 0.25 * np.ones_like(zpos)
        dy= 0.5 * np.ones_like(zpos)
        dz=dat1.values.ravel()
    
        # format dat2 for plotting
        if vmin == None:
            vmin=0
        if vmax == None:
            vmax = np.percentile(dat2,95)
            
        norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
        cmap = matplotlib.cm.get_cmap('jet')
        colors = cmap(norm(dat2.values.ravel()))

        # create figure
        fig = plt.figure(figsize=[14,8])
        ax1 = fig.add_subplot(111, projection='3d')
        
        # make plot
        ax1.grid(True)
        ax1.bar3d(xposM.ravel(),yposM.ravel(),zpos,dx,dy,dz,color=colors,alpha=0.5,edgecolor="None")
        ax1.xaxis.pane.fill = False
        ax1.yaxis.pane.fill = False
        ax1.zaxis.pane.fill = False
        ax1.set_xticks(np.arange(dat1.shape[0])+0.5)
        ax1.set_xticklabels(["%i"%(fw) if (i*100)%200 == 0 else "" for i,fw in zip(range(dat1.shape[0]),flowwidths_cumulative)],rotation=90)  
        ax1.set_yticks(np.arange(dat1.shape[1])+0.5)
        
        # create y-axis ticks
        yticklabels =    pd.date_range(tstart,tend,freq="H")
        if timescale=="hourly":
            ax1.set_yticklabels(["","",""]+[t.strftime("%d %b %y %Hh") if t.hour%tintverval == 0 else "" for t in yticklabels])
        elif timescale=="daily":
            ax1.set_yticklabels(["","",""]+[t.strftime("%d %b %y") if t.day%tintverval == 0 else "" for t in yticklabels])
        elif timescale=="monthly":
            ax1.set_yticklabels(["","",""]+[t.strftime("%b %y") if t.month%tintverval == 0 else "" for t in yticklabels])
        elif timescale=="weekly":
            ax1.set_yticklabels(["","",""]+[t.strftime("%b %y") if t.week%tintverval == 0 else "" for t in yticklabels])

        # set labels of axis
        ax1.set_xlabel("River head --> Outlet [m]", labelpad=20,fontweight='bold')
        ax1.set_ylabel("Time", labelpad=20,fontweight='bold')
        ax1.set_zlabel(lab1,fontweight='bold')        
        ax1.set_title("%s --> %s"%(reach_start,reach_end))
        
        # create colorbar
        ax = fig.add_axes([0.01,0.2,0.1,0.5]) # x,y, lenght, height
        ax.axis("off")
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='20%', pad=0.05)
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm._A = []
        sm
        cb=fig.colorbar(sm,cax=cax, orientation='vertical',extend="max")
        cb.set_label(lab2,fontweight="bold")
        cb.set_clim(vmin=vmin, vmax=vmax)
    
        # save
        
        fname = os.path.join(self.fpath,"branch_efate_and_hydrology_"+fname+".png")
        fig.savefig(fname,dpi=300)      
        
        del fig,dat1,dat2


    def reach_hydrology_and_efate(self,reach,tstart=None,tend=None,
                                  ymax=[0.15,4,0.5,20],maxconc=0.,maxload=0.):
        """
        Plot hydrological and efate state variables of a single reach.
        
        :param reach: name of reach.
        :type reach: str
        :param tstart: Start time.
        :type tstart: datetime
        :param tend: End time.
        :type tend: datetime 
            
        :returns: -
        :rtype: -
        """

        # get simulated data
        res =  self.reach_getTimeseries(reach,variable=None,
                                        tstart=tstart,
                                        tend=tend)
        
       # some plot params
        params=["flow","volume","depth","area"]
        ylabels=["Discharge [m$^3$ sec$^{-1}$]","Volume [m$^3$]","Depth [m]","Wet area [m$^{3}$]"]
        names=["Discharge","Volume","Depth","Wet area"]
        styles=["--","-",":","-."]

        #plot mass balance
        fig = plt.figure(figsize=[12,6])
        host = host_subplot(111, axes_class=AA.Axes)
        plt.subplots_adjust(right=0.75)
       
        # plot subplots right
        offset=0
        for p,l,n,s,ym in zip(params,ylabels,names,styles,ymax):
            par = host.twinx()
            new_fixed_axis = par.get_grid_helper().new_fixed_axis
            par.axis["right"] = new_fixed_axis(loc="right", axes=par,
                                                    offset=(offset, 0))
            par.plot(res.index, res[p], label=l,color="k",linestyle=s)
            offset+=60
            par.set_ylim(0,ym)
            par.set_ylabel(l)
            
        # plot mass flux
        p1, = host.plot(res.index,res["MASS_SW"],color="b",linewidth=3,alpha=.5)
        p1, = host.plot(res.index,res["MASS_SED"],color="orange",linewidth=3,alpha=.5)
        p1, = host.plot(res.index,res["MASS_SED_DEEP"],color="red",linewidth=3,alpha=.5)        
        host.set_xlabel("Time")
        host.set_ylabel("Load [mg]") 
        host.set_ylim(0,maxload)
        #make legend
        plt.legend(loc=0)
        host.axis["bottom"].major_ticklabels.set_rotation(30)
        host.axis["bottom"].major_ticklabels.set_ha("right")
        host.axis["bottom"].label.set_pad(30)
        #save figure
        fname = os.path.join(self.fpath,"reach_hydrology_and_efate_load_"+reach+".png")
        fig.savefig(fname,dpi=300)         
    
        #plot concentration
        fig = plt.figure(figsize=[12,6])
        host = host_subplot(111, axes_class=AA.Axes)
        plt.subplots_adjust(right=0.75)
       
        # plot subplots right
        offset=0
        for p,l,n,s,ym in zip(params,ylabels,names,styles,ymax):
            par = host.twinx()
            new_fixed_axis = par.get_grid_helper().new_fixed_axis
            par.axis["right"] = new_fixed_axis(loc="right", axes=par,
                                                    offset=(offset, 0))
            par.plot(res.index, res[p], label=l,color="k",linestyle=s)
            offset+=60
            par.set_ylim(0,ym)
            par.set_ylabel(l)
    
        # plot mass flux
        p1, = host.plot(res.index,res["PEC_SW"],color="b",linewidth=3,alpha=.5)
        p1, = host.plot(res.index,res["PEC_SED"],color="orange",linewidth=3,alpha=.5)    
        host.set_ylim(0,maxconc)
        host.set_xlabel("Time")
        host.set_ylabel("PEC [$\mu$g L$^{-1}$]")
        host.set_title(reach)
        #make legend
        plt.legend(loc=0)
        host.axis["bottom"].major_ticklabels.set_rotation(30)
        host.axis["bottom"].major_ticklabels.set_ha("right")
        host.axis["bottom"].label.set_pad(30)
        # save figure
        fname = os.path.join(self.fpath,"reach_hydrology_and_efate_conc_"+reach+".png")
        fig.savefig(fname,dpi=300)       
  
    def get_flowwidths(self):
        """
        Get flowwidth of each reach.
        """
        
        # get outlet information
        outlet = self.catchmentlist[self.catchmentlist["component"]=="Outlet"]
        
    
        
        # calcualte flow widths
        flowwidths = []
        for r in self.reaches.index:  
            reach = self.reaches.loc[r]
            # get coords
            x1 = reach.x
            y1 = reach.y
            z1 = reach.z
            if reach.downstream == "Outlet":
                x2 = outlet.x.values[0]
                y2= outlet.y.values[0]
                z2 = outlet.z.values[0]
            else:
                x2 = self.reaches.loc[reach.downstream].x
                y2 = self.reaches.loc[reach.downstream].y
                z2 = self.reaches.loc[reach.downstream].z
                #calculate flow width               
            flowwidths.append((r,np.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)))
        
        # convert to dataframe
        flowwidths = pd.DataFrame(flowwidths,columns=["key","flowwith"])
        flowwidths.set_index("key",inplace=True)

        return flowwidths
    
    def get_slopes(self):
        """
        Get flowwidth of each reach.
        """
        
        # get outlet information
        outlet = self.catchmentlist[self.catchmentlist["component"]=="Outlet"]
        

        # calcualte flow widths
        slopes = []
        for r in self.reaches.index:  
            reach = self.reaches.loc[r]
            # get coords
            x1 = reach.x
            y1 = reach.y
            z1 = reach.z
            if reach.downstream == "Outlet":
                x2 = outlet.x.values[0]
                y2= outlet.y.values[0]
                z2 = outlet.z.values[0]
            else:
                x2 = self.reaches.loc[reach.downstream].x
                y2 = self.reaches.loc[reach.downstream].y
                z2 = self.reaches.loc[reach.downstream].z
                #calculate flow width  
            flowwidth = np.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)
            slopes.append((r,(z2-z1)/(flowwidth)))
        
        # convert to dataframe
        slopes = pd.DataFrame(slopes,columns=["key","slope"])
        slopes.set_index("key",inplace=True)

        return slopes    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    