.. model documentation master file, created by
   sphinx-quickstart on Tue May 29 16:01:42 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive. 

.. figure:: /_static/logo.png
	:width: 7.5cm
	:height: 3cm
	:figclass: align-right

**StepsRiverNetwork**
******************************


Background
==================

StepsRiverNetwork simulates in-stream environmental fate processes of pesticides for
an entire river network of a catchment. Efate processes are calculated for each single 
reach and transport across the entire catchment in an explicit timestep of minutes.

.. figure:: /_static/concept.png
    :align: center
    :figclass: align-center 

For simulating environmental processes the
STEPS1-2-3-4 appraoch by M. Klein (http://publica.fraunhofer.de/documents/N-73445.html) has been implemented. The key features of STEPS-1-2-3-4 are as follow:

* A substance is defined by KOC and DT50.
* The system consists of one water layer and two sediment layer.
* Degradation is calculated by by DT50 in water and sediment.
* Sorption / desorption are calculated by the gradient between water layer and pore water of sediment.

.. figure:: /_static/concept_steps1234.png
    :align: center
    :figclass: align-center 


STEPS1-2-3-4 solely simulates efate processes and hydrological data and compound inputs must be calculated in advance. Having all data at hand efat and solute flux across the whole river network is calculated for each timestep in relation to the followign scheme.

.. figure:: /_static/concept_steps1234_processes.png
    :align: center
    :figclass: align-center 



Input data
==================

All input data files must be located in the respective project folder, except the Runlist.


RunList
-------------------------------

.. figure:: /_static/data_input_RunList.png
    :align: center
    :figclass: align-center 

ReachList
-------------------------------

.. figure:: /_static/data_input_ReachList.png
    :align: center
    :figclass: align-center 

CatchmentList
-------------------------------

.. figure:: /_static/data_input_CatchmentList.png
    :align: center
    :figclass: align-center 

SprayDriftList
-------------------------------

.. figure:: /_static/data_input_SprayDriftList.png
    :align: center
    :figclass: align-center

DrainageList
-------------------------------

.. figure:: /_static/data_input_DrainageList.png
    :align: center
    :figclass: align-center 


HydroList
-------------------------------

The HydroList can be provided by two alternative file formats.

1) CSV

.. figure:: /_static/data_input_reach.png
    :align: center
    :figclass: align-center 

2) HDF5

In this case, the file must called "HydroList.h5". It must contain
three datasets for 'area' (wet surface, m2),
'flow' (m3/day) and 'volume' (m3) according to the parameter described above. Each dataset
must be  a two-dimensional array (shape: [ntime,nreaches]).



Output data
==================

Environmental fate
-------------------------------

1) CSV

A summary CSV-file with a hourly timestep is optionally created which includes the
hydrological and efate data. The hourly summary can be calculated as minimum, maximum or
mean per hour.

.. figure:: /_static/data_output_reach.png
    :align: center
    :figclass: align-center 


2) HDF5

Up to five state variables (MASS_SW,MASS_SED,MASS_SED_DEEP,PEC_SW,PEC_SED)
are saved in a seprated hdf5-file. Each file contains a dataset named by the respective 
parameter (eg. steps_MASS_SW.h5 --> dset="MASS_SW") with a two-dimensional array 
(shape: [ntime,nreaches]).					   

Quick user guide
==================

The model can be executed by using a Python script and calling the related functions:

.. code-block:: python
  :linenos:

  import os
  import argparse
  from pathlib import Path
  from RunFactory import RunFactory

  def get_fdir(subfolder):
    """
    Return current working directoy and add subfolder to path.
    """
    fdir = os.path.abspath(os.path.join(
                    os.path.dirname(Path(__file__).parent),
                    *subfolder))
    return fdir

  FLAGS = None

  if __name__ == "__main__":
    
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




... or by using a batchfile by defining the variables '--folder' and '--runlist' which 
area the the project folder and the name of the project:

.. code-block:: bash
  :linenos:

  @echo off
  set script=%cd%/bin/main.py
  set python=%cd%/bin/Python/python.exe
  call %python% %script% --folder c:/projects/test/ --runlist testrun
  pause


About
==================

The tool is a development in a project by Bayer AG and knoell Germany GmbH.


Sebastian Multsch :sup:`1` , Stefan Reichenberger :sup:`1` ,Florian Krebs :sup:`1` , Thorsten Schad :sup:`2` 

:sup:`1` `knoell Germany GmbH <https://www.knoellconsult.com/enf>`_ 

:sup:`2` `Bayer AG, Research & Development, Crop Science <https://www.cropscience.bayer.de/>`_ 
