## Table of Contents
* [About the project](#about-the-project)
  * [Built With](#built-with)
* [Getting Started](#getting-started)
  * [Prerequisites](#prerequisites)
  * [Installation](#installation)
* [Usage](#usage)
  * [Inputs](#inputs)
  * [Outputs](#outputs)
* [Roadmap](#roadmap)
* [Contributing](#contributing)
* [License](#license)
* [Contact](#contact)
* [Acknowledgements](#acknowledgements)


## About the project
The component encapsulating the Steps environmental fate module.  
This is an automatically generated documentation based on the available code and in-line documentation. The current
version of this document is from 2021-08-27.  

### Built with
* Landscape Model core version 1.5.9
* River network version of STEPS1234 version 0.93 (see `module\documentation\html\index.html` for details)


## Getting Started
The component can be used in any Landscape Model based on core version 1.5.9 or newer. See the Landscape Model
core's `README` for general tips on how to add a component to a Landscape Model.

### Prerequisites
A model developer that wants to add the `StepsRiverNetwork` component to a Landscape Model needs to set up the general 
structure for a Landscape Model first. See the Landscape Model core's `README` for details on how to do so.

### Installation
1. Copy the `StepsRiverNetwork` component into the `model\variant` sub-folder.
2. Make use of the component by including it into the model composition using `module=StepsRiverNetwork` and 
   `class=StepsRiverNetwork`. 


## Usage
The following gives a sample configuration of the `StepsRiverNetwork` component. See [inputs](#inputs) and 
[outputs](#outputs) for further details on the component's interface.
```xml
<StepsRiverNetwork module="StepsRiverNetwork" class="StepsRiverNetwork" enabled="$(RunStepsRiverNetwork)">
<ProcessingPath>$(_MCS_BASE_DIR_)\$(_MC_NAME_)\processing\fate\steps</ProcessingPath>
    <Hydrography>
<FromOutput component="LandscapeScenario" output="hydrography" />
    </Hydrography>
<Catchment>$(:Catchment)</Catchment>
    <WaterDischarge>
        <FromOutput component="Hydrology" output="Flow" />
</WaterDischarge>
    <TimeSeriesStart>
        <FromOutput component="Hydrology" output="TimeSeriesStart" />
</TimeSeriesStart>
    <ReachesHydrology>
        <FromOutput component="Hydrology" output="Reaches" />
</ReachesHydrology>
    <WaterVolume>
        <FromOutput component="Hydrology" output="Volume" />
    </WaterVolume>
<WetSurfaceArea>
        <FromOutput component="Hydrology" output="Area" />
    </WetSurfaceArea>
    <DriftDeposition>
<FromOutput component="DepositionToReach" output="Deposition" />
    </DriftDeposition>
    <ReachesDrift>
<FromOutput component="DepositionToReach" output="Reaches" />
    </ReachesDrift>
    <MolarMass type="float"
unit="g/mol">$(MolarMass)</MolarMass>
    <DT50sw type="float" unit="d">$(DT50sw)</DT50sw>
    <DT50sed type="float"
unit="d">$(DT50sed)</DT50sed>
    <KOC type="float" unit="l/kg">$(KOC)</KOC>
    <Temp0 type="float"
unit="&#176;C">$(Temp0)</Temp0>
    <Q10 type="float" unit="1">$(Q10)</Q10>
    <PlantUptake type="float"
unit="1">$(PlantUptake)</PlantUptake>
    <QFac type="float" unit="1">$(QFac)</QFac>
    <ThresholdSW type="float"
unit="mg/m&#179;">$(ThresholdSW)</ThresholdSW>
    <ThresholdSediment type="float"
unit="mg/kg">$(ThresholdSediment)</ThresholdSediment>
</StepsRiverNetwork>
```

### Inputs
#### ProcessingPath
The working directory for the module. It is used for all files prepared as module inputs
or generated as (temporary) module outputs.  
`ProcessingPath` expects its values to be of type `str`.
Values of the `ProcessingPath` input may not have a physical unit.
Values have to refer to the `global` scale.

#### Hydrography
The spatial delineation of the hydrographic features in the simulated landscape. This
input basically represents the flow-lines used during preparation of the hydrology. The hydrography is
consistently for all components of the Landscape Model subdivided into individual segments (*reaches*).  
`Hydrography` expects its values to be of type `str`.
Values of the `Hydrography` input may not have a physical unit.
Values have to refer to the `global` scale.

#### Catchment
A file path to a CSV file detailing the hydrographic properties of the entire catchment
depicted by hydrographic the scenario. This file is usually provided by the scenario developer (if
usage of StepsRiverNetwork is supported by the scenario) and is made available as a project macro.  
`Catchment` expects its values to be of type `str`.
Values of the `Catchment` input may not have a physical unit.
Values have to refer to the `global` scale.

#### WaterDischarge
The entire water discharge of this reach into the next downstream reach.  
`WaterDischarge` expects its values to be of type `ndarray`.
The physical unit of the `WaterDischarge` input values is `m³/d`.
Values have to refer to the `time/hour, space/reach` scale.

#### TimeSeriesStart
The first time step for which input data is provided. This is also the time step of where
the StepsRiverNetwork simulation starts.  
`TimeSeriesStart` expects its values to be of type `datetime`.
Values of the `TimeSeriesStart` input may not have a physical unit.
Values have to refer to the `global` scale.

#### ReachesHydrology
The numeric identifiers for individual reaches (in the order of the hydrological inputs).  
`ReachesHydrology` expects its values to be of type `ndarray`.
Values of the `ReachesHydrology` input may not have a physical unit.
Values have to refer to the `space/reach` scale.

#### WaterVolume
The amount of water contained by a reach.  
`WaterVolume` expects its values to be of type `ndarray`.
The physical unit of the `WaterVolume` input values is `m³`.
Values have to refer to the `time/hour, space/reach` scale.

#### WetSurfaceArea
The surface area of a reach.  
`WetSurfaceArea` expects its values to be of type `ndarray`.
The physical unit of the `WetSurfaceArea` input values is `m²`.
Values have to refer to the `time/hour, space/reach` scale.

#### DriftDeposition
The average drift deposition onto the surface of a water body.  
`DriftDeposition` expects its values to be of type `ndarray`.
The physical unit of the `DriftDeposition` input values is `mg/m²`.
Values have to refer to the `time/day, space/reach` scale.

#### ReachesDrift
The numeric identifiers for individual reaches (in the order of the `DriftDeposition` 
input) that apply scenario-wide.  
`ReachesDrift` expects its values to be of type `ndarray`.
Values of the `ReachesDrift` input may not have a physical unit.
Values have to refer to the `space/reach` scale.

#### MolarMass
The molar mass of the substance depositing at the water body surface.  
`MolarMass` expects its values to be of type `float`.
The physical unit of the `MolarMass` input values is `g/mol`.

#### DT50sw
The half-life transformation time in water of the substance depositing at the water body 
surface.  
`DT50sw` expects its values to be of type `float`.
The physical unit of the `DT50sw` input values is `d`.

#### DT50sed
The half-life transformation time in sediment of the substance depositing at the water 
body surface.  
`DT50sed` expects its values to be of type `float`.
The physical unit of the `DT50sed` input values is `d`.

#### KOC
The coefficient for equilibrium adsorption in sediment of the substance depositing at 
the water body surface.  
`KOC` expects its values to be of type `float`.
The physical unit of the `KOC` input values is `l/kg`.

#### Temp0
The reference temperature to which the physical and chemical substance values apply.  
`Temp0` expects its values to be of type `float`.
The physical unit of the `Temp0` input values is `°C`.

#### Q10
The temperature coefficient for chemical reactions of the deposited substance.  
`Q10` expects its values to be of type `float`.
The physical unit of the `Q10` input values is `1`.

#### PlantUptake
The fraction of pesticide that is taken up by plants.  
`PlantUptake` expects its values to be of type `float`.
The physical unit of the `PlantUptake` input values is `1`.

#### QFac
The QFac parameter is not documented in the module documentation.  
`QFac` expects its values to be of type `float`.
The physical unit of the `QFac` input values is `1`.

#### ThresholdSW
The minimum surface water concentration that is reported.  
`ThresholdSW` expects its values to be of type `float`.
The physical unit of the `ThresholdSW` input values is `mg/m³`.

#### ThresholdSediment
The minimum sediment concentration that is reported.  
`ThresholdSediment` expects its values to be of type `float`.
The physical unit of the `ThresholdSediment` input values is `mg/kg`.

### Outputs
#### PEC_SW
The modelled concentration in the water phase.  
Values are expectedly of type `ndarray`.
Value representation is in a 2-dimensional array.
Dimension 1 spans the number of simulated hours.
Dimension 2 spans the number of simulated reaches.
Chunking of the array is for fast retrieval of time series.
Values expectedly have a unit of `mg/m³`.
Individual array elements have a type of `float`.
The values apply to the following scale: `time/hour, space/base_geometry`.
#### MASS_SW
The modelled substance mass in the water phase.  
Values are expectedly of type `ndarray`.
Value representation is in a 2-dimensional array.
Dimension 1 spans the number of simulated hours.
Dimension 2 spans the number of simulated reaches.
Chunking of the array is for fast retrieval of time series.
Values expectedly have a unit of `mg`.
Individual array elements have a type of `float`.
The values apply to the following scale: `time/hour, space/base_geometry`.
#### MASS_SED
The modelled substance mass in sediment.  
Values are expectedly of type `ndarray`.
Value representation is in a 2-dimensional array.
Dimension 1 spans the number of simulated hours.
Dimension 2 spans the number of simulated reaches.
Chunking of the array is for fast retrieval of time series.
Values expectedly have a unit of `mg`.
Individual array elements have a type of `float`.
The values apply to the following scale: `time/hour, space/base_geometry`.
#### MASS_SED_DEEP
The modelled substance mass in deep sediment.  
Values are expectedly of type `ndarray`.
Value representation is in a 2-dimensional array.
Dimension 1 spans the number of simulated hours.
Dimension 2 spans the number of simulated reaches.
Chunking of the array is for fast retrieval of time series.
Values expectedly have a unit of `mg`.
Individual array elements have a type of `float`.
The values apply to the following scale: `time/hour, space/base_geometry`.
#### PEC_SED
The modelled concentration in sediment.  
Values are expectedly of type `ndarray`.
Value representation is in a 2-dimensional array.
Dimension 1 spans the number of simulated hours.
Dimension 2 spans the number of simulated reaches.
Chunking of the array is for fast retrieval of time series.
Values expectedly have a unit of `mg/m³`.
Individual array elements have a type of `float`.
The values apply to the following scale: `time/hour, space/base_geometry`.
#### Reaches
The numerical identifiers of the reaches in the order of the other outputs.  
Values are expectedly of type `list[int]`.
The values apply to the following scale: `space/reach`.


## Roadmap
The following changes will be part of future `StepsRiverNetwork` versions:
* z-value precision ([#1](https://gitlab.bayer.com/aqrisk-landscape/stepsrivernetwork-component/-/issues/1))


## Contributing
Contributions are welcome. Please contact the authors (see [Contact](#contact)). Also consult the `CONTRIBUTING` 
document for more information.


## License
Distributed under the CC0 License. See `LICENSE` for more information.


## Contact
Sascha Bub (component) - sascha.bub@gmx.de  
Thorsten Schad (component) - thorsten.schad@bayer.com  
Sebastian Multsch (module) - smultsch@knoell.com  


## Acknowledgements
* [GDAL](https://pypi.org/project/GDAL)  
* [h5py](https://www.h5py.org)  
* [NumPy](https://numpy.org)  
