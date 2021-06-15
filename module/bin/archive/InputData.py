# -*- coding: utf-8 -*-
"""
Created on Mon Aug  6 15:25:21 2018

@author: smh
"""

import os
from xml.dom import minidom
from Parameter import ParameterList

class InputData(object):
    """
    Iterable class which holds a set of data records
    """
    def __init__(self,fpath,sep=","):
        
        # read input data information from xml
        key = [i.attributes['name'].value for i in self.__read_info(
                os.sep.join(os.path.abspath(__file__).split("\\")[:-1]),"input_tables.xml")]
        value = [ParameterList(fpath,fname,sep) for fname in key]
         
        for (key, value) in zip(key, value):
           self.__dict__[key] = value
    
    def __setattr__(self, name, value):
        raise Exception("Read only")

    def __read_info(self,fpath,fname):
        # parse an xml file by name
        mydoc = minidom.parse(os.path.join(fpath,fname))
        items = mydoc.getElementsByTagName('item')
        return items
    
    
if __name__ == "__main__":

    fpath = "c:/_STEPSRivernetwork/projects/Rummen_subCatch_20reaches/hydro_v01_medium/"
    
