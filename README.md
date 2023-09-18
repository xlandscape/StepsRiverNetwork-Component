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

The component encapsulating the StepsRiverNetwork environmental fate module. StepsRiverNetwork simulates in-stream
environmental fate processes of pesticides for an entire river network of a catchment. Environmental fate processes
are calculated for each single reach, and transport across the entire catchment is reported in an explicit timestep
of one hour.  
This is an automatically generated documentation based on the available code and in-line documentation. The current
version of this document is from 2023-09-18.

### Built with

* Landscape Model core version 1.15.5
* River network version of STEPS1234 version 0.93 (see `module\documentation\html\index.html` for details)

## Getting Started

The component can be used in any Landscape Model based on core version 1.15.5 or newer. See the Landscape
Model core's `README` for general tips on how to add a component to a Landscape Model.

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
<ProcessingPath scales="global">$(_MCS_BASE_DIR_)\$(_MC_NAME_)\processing\fate\steps</ProcessingPath>
    <Catchment
scales="global">$(:Catchment)</Catchment>
    <WaterDischarge>
        <FromOutput component="Hydrology" output="Flow"
/>
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
    <MolarMass type="float"
unit="g/mol" scales="global">$(MolarMass)</MolarMass>
    <DT50sw type="float" unit="d"
scales="global">$(DT50sw)</DT50sw>
    <DT50sed type="float" unit="d" scales="global">$(DT50sed)</DT50sed>
    <KOC
type="float" unit="l/kg" scales="global">$(KOC)</KOC>
    <Temp0 type="float" unit="&#176;C"
scales="global">$(Temp0)</Temp0>
    <Q10 type="float" unit="1" scales="global">$(Q10)</Q10>
    <PlantUptake
type="float" unit="1" scales="global">$(PlantUptake)</PlantUptake>
    <QFac type="float" unit="1"
scales="global">$(QFac)</QFac>
    <ThresholdSW type="float" unit="mg/m&#179;"
scales="global">$(ThresholdSW)</ThresholdSW>
    <ThresholdSediment type="float" unit="mg/kg"
scales="global">$(ThresholdSediment)</ThresholdSediment>
    <HydrographyGeometries>
        <FromOutput
component="LandscapeScenario" output="hydrography_geom" />
    </HydrographyGeometries>
    <DownstreamReach>
<FromOutput component="LandscapeScenario" output="hydrography_downstream" />
    </DownstreamReach>
    <InitialDepth>
<FromOutput component="LandscapeScenario" output="hydrography_initial_depth" />
    </InitialDepth>
    <Manning>
<FromOutput component="LandscapeScenario" output="hydrography_manning" />
    </Manning>
    <BankSlope>
<FromOutput component="LandscapeScenario" output="hydrography_bank_slope" />
    </BankSlope>
    <Width>
<FromOutput component="LandscapeScenario" output="hydrography_width" />
    </Width>
    <Shape>
        <FromOutput
component="LandscapeScenario" output="hydrography_shape" />
    </Shape>
    <BulkDensity>
        <FromOutput
component="LandscapeScenario" output="hydrography_bulk_density" />
    </BulkDensity>
    <Porosity>
        <FromOutput
component="LandscapeScenario" output="hydrography_porosity" />
    </Porosity>
    <OrganicContent>
        <FromOutput
component="LandscapeScenario" output="hydrography_organic_content" />
    </OrganicContent>
    <SedimentDepth1stLayer>
<FromOutput component="LandscapeScenario" output="hydrography_sediment_layer_1_depth" />
    </SedimentDepth1stLayer>
<SedimentDepth2ndLayer>
        <FromOutput component="LandscapeScenario" output="hydrography_sediment_layer_2_depth" />
</SedimentDepth2ndLayer>
</StepsRiverNetwork>
```

### Inputs

#### ProcessingPath

The working directory for the module. It is used for all files prepared as module inputs or generated as (temporary)
module outputs.
`ProcessingPath` expects its values to be of type `str`.
Values of the `ProcessingPath` input may not have a physical unit.
Values have to refer to the `global` scale.

#### Catchment

A file path to a CSV file detailing the hydrographic properties of the entire catchment depicted by a hydrographic
scenario. This file is usually provided by the scenario developer (if usage of StepsRiverNetwork is supported by the
scenario) and is made available as a project macro. See the module documentation for details on the format.
`Catchment` expects its values to be of type `str`.
Values of the `Catchment` input may not have a physical unit.
Values have to refer to the `global` scale.

#### WaterDischarge

`WaterDischarge` expects its values to be of type `ndarray`.
The physical unit of the `WaterDischarge` input values is `m³/d`.
Values have to refer to the `time/hour, space/reach` scale.

#### TimeSeriesStart

The first time step for which input data is provided. This is also the time step of where the StepsRiverNetwork
simulation starts. This input will be removed in a future version of the `StepsRiverNetwork` component.
`TimeSeriesStart` expects its values to be of type `datetime`.
Values of the `TimeSeriesStart` input may not have a physical unit.
Values have to refer to the `global` scale.

#### ReachesHydrology

The numeric identifiers for individual reaches (in the order of the hydrological inputs). his input will be removed in a
future version of the `StepsRiverNetwork` component.
`ReachesHydrology` expects its values to be of type `ndarray`.
Values of the `ReachesHydrology` input may not have a physical unit.
Values have to refer to the `space/reach` scale.

#### WaterVolume

`WaterVolume` expects its values to be of type `ndarray`.
The physical unit of the `WaterVolume` input values is `m³`.
Values have to refer to the `time/hour, space/reach` scale.

#### WetSurfaceArea

`WetSurfaceArea` expects its values to be of type `ndarray`.
The physical unit of the `WetSurfaceArea` input values is `m²`.
Values have to refer to the `time/hour, space/reach` scale.

#### DriftDeposition

`DriftDeposition` expects its values to be of type `ndarray`.
The physical unit of the `DriftDeposition` input values is `mg/m²`.
Values have to refer to the `time/day, space/reach` scale.

#### MolarMass

`MolarMass` expects its values to be of type `float`.
The physical unit of the `MolarMass` input values is `g/mol`.
Values have to refer to the `global` scale.

#### DT50sw

`DT50sw` expects its values to be of type `float`.
The physical unit of the `DT50sw` input values is `d`.
Values have to refer to the `global` scale.

#### DT50sed

`DT50sed` expects its values to be of type `float`.
The physical unit of the `DT50sed` input values is `d`.
Values have to refer to the `global` scale.

#### KOC

`KOC` expects its values to be of type `float`.
The physical unit of the `KOC` input values is `l/kg`.
Values have to refer to the `global` scale.

#### Temp0

`Temp0` expects its values to be of type `float`.
The physical unit of the `Temp0` input values is `°C`.
Values have to refer to the `global` scale.

#### Q10

`Q10` expects its values to be of type `float`.
The physical unit of the `Q10` input values is `1`.
Values have to refer to the `global` scale.

#### PlantUptake

`PlantUptake` expects its values to be of type `float`.
The physical unit of the `PlantUptake` input values is `1`.
Values have to refer to the `global` scale.

#### QFac

`QFac` expects its values to be of type `float`.
The physical unit of the `QFac` input values is `1`.
Values have to refer to the `global` scale.

#### ThresholdSW

`ThresholdSW` expects its values to be of type `float`.
The physical unit of the `ThresholdSW` input values is `mg/m³`.
Values have to refer to the `global` scale.

#### ThresholdSediment

`ThresholdSediment` expects its values to be of type `float`.
The physical unit of the `ThresholdSediment` input values is `mg/kg`.
Values have to refer to the `global` scale.

#### HydrographyGeometries

The geometries of individual water body segments (reaches) in WKB representation. This input will be removed in a future
version of the `StepsRiverNetwork` component.
`HydrographyGeometries` expects its values to be of type `list`.
Values of the `HydrographyGeometries` input may not have a physical unit.
Values have to refer to the `space/reach` scale.

#### DownstreamReach

`DownstreamReach` expects its values to be of type `list`.
Values of the `DownstreamReach` input may not have a physical unit.
Values have to refer to the `space/reach` scale.

#### InitialDepth

`InitialDepth` expects its values to be of type `list`.
The physical unit of the `InitialDepth` input values is `m`.
Values have to refer to the `space/reach` scale.

#### Manning

`Manning` expects its values to be of type `list`.
The physical unit of the `Manning` input values is `1`.
Values have to refer to the `space/reach` scale.

#### BankSlope

`BankSlope` expects its values to be of type `list`.
The physical unit of the `BankSlope` input values is `1`.
Values have to refer to the `space/reach` scale.

#### Width

`Width` expects its values to be of type `list`.
The physical unit of the `Width` input values is `m`.
Values have to refer to the `space/reach` scale.

#### Shape

`Shape` expects its values to be of type `list`.
Values of the `Shape` input may not have a physical unit.
Values have to refer to the `space/reach` scale.
Allowed values are: `TriangularReach`, `RectangularReach`, `SWATReachType`.

#### BulkDensity

`BulkDensity` expects its values to be of type `list`.
The physical unit of the `BulkDensity` input values is `kg/m³`.
Values have to refer to the `space/reach` scale.

#### Porosity

`Porosity` expects its values to be of type `list`.
The physical unit of the `Porosity` input values is `m³/m³`.
Values have to refer to the `space/reach` scale.

#### OrganicContent

`OrganicContent` expects its values to be of type `list`.
The physical unit of the `OrganicContent` input values is `g/g`.
Values have to refer to the `space/reach` scale.

#### SedimentDepth1stLayer

`SedimentDepth1stLayer` expects its values to be of type `list`.
The physical unit of the `SedimentDepth1stLayer` input values is `m`.
Values have to refer to the `space/reach` scale.

#### SedimentDepth2ndLayer

`SedimentDepth2ndLayer` expects its values to be of type `list`.
The physical unit of the `SedimentDepth2ndLayer` input values is `m`.
Values have to refer to the `space/reach` scale.

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
The values apply to the following scale: `time/hour, space/reach`.
#### MASS_SW
The modelled substance mass in the water phase.  
Values are expectedly of type `ndarray`.
Value representation is in a 2-dimensional array.
Dimension 1 spans the number of simulated hours.
Dimension 2 spans the number of simulated reaches.
Chunking of the array is for fast retrieval of time series.
Values expectedly have a unit of `mg`.
Individual array elements have a type of `float`.
The values apply to the following scale: `time/hour, space/reach`.
#### MASS_SED
The modelled substance mass in sediment.  
Values are expectedly of type `ndarray`.
Value representation is in a 2-dimensional array.
Dimension 1 spans the number of simulated hours.
Dimension 2 spans the number of simulated reaches.
Chunking of the array is for fast retrieval of time series.
Values expectedly have a unit of `mg`.
Individual array elements have a type of `float`.
The values apply to the following scale: `time/hour, space/reach`.
#### MASS_SED_DEEP
The modelled substance mass in deep sediment.  
Values are expectedly of type `ndarray`.
Value representation is in a 2-dimensional array.
Dimension 1 spans the number of simulated hours.
Dimension 2 spans the number of simulated reaches.
Chunking of the array is for fast retrieval of time series.
Values expectedly have a unit of `mg`.
Individual array elements have a type of `float`.
The values apply to the following scale: `time/hour, space/reach`.
#### PEC_SED
The modelled concentration in sediment.  
Values are expectedly of type `ndarray`.
Value representation is in a 2-dimensional array.
Dimension 1 spans the number of simulated hours.
Dimension 2 spans the number of simulated reaches.
Chunking of the array is for fast retrieval of time series.
Values expectedly have a unit of `mg/m³`.
Individual array elements have a type of `float`.
The values apply to the following scale: `time/hour, space/reach`.
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
