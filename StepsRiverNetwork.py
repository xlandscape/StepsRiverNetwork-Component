"""
Component for the Steps environmental fate module.
"""
import datetime
import h5py
import numpy as np
from osgeo import ogr
import os
import shutil
import base
import attrib


class StepsRiverNetwork(base.Component):
    """
    The component encapsulating the Steps environmental fate module.
    """
    # RELEASES
    VERSION = base.VersionCollection(
        base.VersionInfo("2.1.2", "2021-12-10"),
        base.VersionInfo("2.1.1", "2021-11-18"),
        base.VersionInfo("2.1.0", "2021-10-20"),
        base.VersionInfo("2.0.7", "2021-10-19"),
        base.VersionInfo("2.0.6", "2021-10-12"),
        base.VersionInfo("2.0.5", "2021-10-11"),
        base.VersionInfo("2.0.4", "2021-09-01"),
        base.VersionInfo("2.0.3", "2021-08-27"),
        base.VersionInfo("2.0.2", "2021-07-19"),
        base.VersionInfo("2.0.1", "2020-12-03"),
        base.VersionInfo("2.0.0", "2020-10-22"),
        base.VersionInfo("1.3.35", "2020-08-12"),
        base.VersionInfo("1.3.33", "2020-07-30"),
        base.VersionInfo("1.3.29", "2020-06-15"),
        base.VersionInfo("1.3.27", "2020-05-20"),
        base.VersionInfo("1.3.26", "2020-04-09"),
        base.VersionInfo("1.3.24", "2020-04-02"),
        base.VersionInfo("1.3.16", "2020-02-11"),
        base.VersionInfo("1.3.10", "2020-01-23"),
        base.VersionInfo("1.3.7", "2020-01-10"),
        base.VersionInfo("1.3.3", "2019-12-15"),
        base.VersionInfo("1.2.38", None),
        base.VersionInfo("1.2.37", None),
        base.VersionInfo("1.2.36", None),
        base.VersionInfo("1.2.4", None),
        base.VersionInfo("1.2.3", None),
        base.VersionInfo("1.1.1", None)
    )

    # AUTHORS
    VERSION.authors.extend((
        "Sascha Bub (component) - sascha.bub@gmx.de",
        "Thorsten Schad (component) - thorsten.schad@bayer.com",
        "Sebastian Multsch (module) - smultsch@knoell.com"
    ))

    # ACKNOWLEDGEMENTS
    VERSION.acknowledgements.extend((
        "[GDAL](https://pypi.org/project/GDAL)",
        "[h5py](https://www.h5py.org)",
        "[NumPy](https://numpy.org)"
    ))

    # ROADMAP
    VERSION.roadmap.extend((
        "z-value precision ([#1](https://gitlab.bayer.com/aqrisk-landscape/stepsrivernetwork-component/-/issues/1))",
    ))

    # CHANGELOG
    # noinspection SpellCheckingInspection
    VERSION.added("1.2.3", "`components.StepsRivernetwork` component")
    VERSION.added("1.2.4", "`components.CmfHydrology` hydrology Python library only loaded when needed")
    VERSION.changed("1.2.36", "`components.CmfHydrology` replaced by `components.StepsRiverNetwork` ")
    # noinspection SpellCheckingInspection
    VERSION.changed("1.2.37", "`components.StepsRivernetwork` module updated to version 0.9.2")
    # noinspection SpellCheckingInspection
    VERSION.changed("1.2.37", "`components.StepsRivernetwork` thresholds as inputs")
    # noinspection SpellCheckingInspection
    VERSION.fixed("1.3.3", "increased numeric precision in `components.StepsRivernetwork` (preliminary)")
    # noinspection SpellCheckingInspection
    VERSION.changed("1.3.7", "`components.StepsRivernetwork` module updated to version 0.9.3")
    # noinspection SpellCheckingInspection
    VERSION.fixed("1.3.10", "Spatial referencing in `components.StepsRivernetwork`")
    # noinspection SpellCheckingInspection
    VERSION.changed("1.3.16", "Substance parameterization in `components.StepsRivernetwork` changed")
    # noinspection SpellCheckingInspection
    VERSION.changed("1.3.24", "`components.StepsRivernetwork` uses base function to call module")
    # noinspection SpellCheckingInspection
    VERSION.fixed("1.3.26", "`components.StepsRivernetwork` reach order")
    # noinspection SpellCheckingInspection
    VERSION.changed("1.3.27", "`components.StepsRivernetwork` specifies scales")
    # noinspection SpellCheckingInspection
    VERSION.fixed("1.3.29", "Input slicing in `components.StepsRivernetwork`")
    # noinspection SpellCheckingInspection
    VERSION.changed("1.3.33", "`components.StepsRivernetwork` checks input types strictly")
    # noinspection SpellCheckingInspection
    VERSION.changed("1.3.33", "`components.StepsRivernetwork` checks for physical units")
    # noinspection SpellCheckingInspection
    VERSION.changed("1.3.33", "`components.StepsRivernetwork` reports physical units to the data store")
    # noinspection SpellCheckingInspection
    VERSION.changed("1.3.33", "`components.StepsRivernetwork` checks for scales")
    # noinspection SpellCheckingInspection
    VERSION.changed(
        "1.3.35", "`components.StepsRivernetwork` receives processing path as home path environment variable")
    VERSION.changed("2.0.0", "First independent release")
    VERSION.added("2.0.1", "Changelog and release history")
    VERSION.changed("2.0.2", "Spellings")
    VERSION.changed("2.0.2", "Changelog uses markdown")
    VERSION.added("2.0.3", "Base documentation")
    VERSION.changed("2.0.4", "ogr module import")
    VERSION.changed("2.0.4", "Acknowledged default access mode for HDF files")
    VERSION.changed("2.0.5", "Replaced legacy format strings by f-strings")
    VERSION.changed("2.0.6", "Switched to Google docstring style")
    VERSION.changed("2.0.7", "Specified working directory for module")
    VERSION.changed("2.1.0", "Replaced shapefile input")
    VERSION.changed("2.1.1", "Removed reaches inputs")
    VERSION.changed("2.1.1", "Reports element names of outputs")
    VERSION.changed("2.1.2", "Specifies offset of outputs")

    def __init__(self, name, observer, store):
        super(StepsRiverNetwork, self).__init__(name, observer, store)
        self._module = base.Module(
            "River network version of STEPS1234", "0.93", r"module\documentation\html\index.html")
        # noinspection SpellCheckingInspection
        self._inputs = base.InputContainer(self, [
            base.Input(
                "ProcessingPath",
                (attrib.Class(str, 1), attrib.Unit(None, 1), attrib.Scales("global", 1)),
                self.default_observer,
                description="""The working directory for the module. It is used for all files prepared as module inputs
                or generated as (temporary) module outputs."""
            ),
            base.Input(
                "Hydrography",
                (attrib.Class(str, 1), attrib.Unit(None, 1), attrib.Scales("global", 1)),
                self.default_observer,
                description="""The spatial delineation of the hydrographic features in the simulated landscape. This
                consistently for all components of the Landscape Model subdivided into individual segments (*reaches*).
                """
            ),
            base.Input(
                "Catchment",
                (attrib.Class(str, 1), attrib.Unit(None, 1), attrib.Scales("global", 1)),
                self.default_observer,
                description="""A file path to a CSV file detailing the hydrographic properties of the entire catchment
                depicted by hydrographic the scenario. This file is usually provided by the scenario developer (if
                usage of StepsRiverNetwork is supported by the scenario) and is made available as a project macro."""
            ),
            base.Input(
                "WaterDischarge",
                (
                    attrib.Class(np.ndarray, 1),
                    attrib.Unit("m³/d", 1),
                    attrib.Scales("time/hour, space/reach", 1)
                ),
                self.default_observer,
                description="The entire water discharge of this reach into the next downstream reach."
            ),
            base.Input(
                "TimeSeriesStart",
                (attrib.Class(datetime.datetime, 1), attrib.Unit(None, 1), attrib.Scales("global", 1)),
                self.default_observer,
                description="""The first time step for which input data is provided. This is also the time step of where
                the StepsRiverNetwork simulation starts."""
            ),
            base.Input(
                "ReachesHydrology",
                (attrib.Class(np.ndarray, 1), attrib.Unit(None, 1), attrib.Scales("space/reach", 1)),
                self.default_observer,
                description="The numeric identifiers for individual reaches (in the order of the hydrological inputs)."
            ),
            base.Input(
                "WaterVolume",
                (
                    attrib.Class(np.ndarray, 1),
                    attrib.Unit("m³", 1),
                    attrib.Scales("time/hour, space/reach", 1)
                ),
                self.default_observer,
                description="The amount of water contained by a reach."
            ),
            base.Input(
                "WetSurfaceArea",
                (
                    attrib.Class(np.ndarray, 1),
                    attrib.Unit("m²", 1),
                    attrib.Scales("time/hour, space/reach", 1)
                ),
                self.default_observer,
                description="The surface area of a reach."
            ),
            base.Input(
                "DriftDeposition",
                (
                    attrib.Class(np.ndarray, 1),
                    attrib.Unit("mg/m²", 1),
                    attrib.Scales("time/day, space/reach", 1)
                ),
                self.default_observer,
                description="The average drift deposition onto the surface of a water body."
            ),
            base.Input(
                "DrainageLoad",
                (
                    attrib.Class(np.ndarray, 1),
                    attrib.Unit("mg/m2/h", 1),
                    attrib.Scales("time/hour, space/base_geometry", 1)
                ),
                self.default_observer
            ),
            base.Input(
                "ReachesDrift",
                (attrib.Class(np.ndarray, 1), attrib.Unit(None, 1), attrib.Scales("space/reach", 1)),
                self.default_observer,
                description="""The numeric identifiers for individual reaches (in the order of the `DriftDeposition` 
                input) that apply scenario-wide."""
            ),
            base.Input(
                "MolarMass",
                (attrib.Class(float, 1), attrib.Unit("g/mol", 1)),
                self.default_observer,
                description="The molar mass of the substance depositing at the water body surface."
            ),
            base.Input(
                "DT50sw",
                (attrib.Class(float, 1), attrib.Unit("d", 1)),
                self.default_observer,
                description="""The half-life transformation time in water of the substance depositing at the water body 
                surface."""
            ),
            base.Input(
                "DT50sed",
                (attrib.Class(float, 1), attrib.Unit("d", 1)),
                self.default_observer,
                description="""The half-life transformation time in sediment of the substance depositing at the water 
                body surface."""
            ),
            base.Input(
                "KOC",
                (attrib.Class(float, 1), attrib.Unit("l/kg", 1)),
                self.default_observer,
                description="""The coefficient for equilibrium adsorption in sediment of the substance depositing at 
                the water body surface."""
            ),
            base.Input(
                "Temp0",
                (attrib.Class(float, 1), attrib.Unit("°C", 1)),
                self.default_observer,
                description="The reference temperature to which the physical and chemical substance values apply."
            ),
            base.Input(
                "Q10",
                (attrib.Class(float, 1), attrib.Unit("1", 1)),
                self.default_observer,
                description="The temperature coefficient for chemical reactions of the deposited substance."
            ),
            base.Input(
                "PlantUptake",
                (attrib.Class(float, 1), attrib.Unit("1", 1)),
                self.default_observer,
                description="The fraction of pesticide that is taken up by plants."
            ),
            base.Input(
                "QFac",
                (attrib.Class(float, 1), attrib.Unit("1", 1)),
                self.default_observer,
                description="The QFac parameter is not documented in the module documentation."
            ),
            base.Input(
                "ThresholdSW",
                (attrib.Class(float, 1), attrib.Unit("mg/m³", 1)),
                self.default_observer,
                description="The minimum surface water concentration that is reported."
            ),
            base.Input(
                "ThresholdSediment",
                (attrib.Class(float, 1), attrib.Unit("mg/kg", 1)),
                self.default_observer,
                description="The minimum sediment concentration that is reported."
            ),
            base.Input(
                "HydrographyGeometries",
                (attrib.Class(list[bytes]), attrib.Unit(None), attrib.Scales("space/base_geometry")),
                self.default_observer,
                description="The geometries of individual water body segments (reaches) in WKB representation."
            ),
            base.Input(
                "DownstreamReach",
                (attrib.Class(list[str]), attrib.Unit(None), attrib.Scales("space/base_geometry")),
                self.default_observer,
                description="The identifier of the reach that is located downstream of the current reach."
            ),
            base.Input(
                "InitialDepth",
                (attrib.Class(list[float]), attrib.Unit("m"), attrib.Scales("space/base_geometry")),
                self.default_observer,
                description="The initial water depth of the current reach."
            ),
            base.Input(
                "Manning",
                (attrib.Class(list[float]), attrib.Unit("1"), attrib.Scales("space/base_geometry")),
                self.default_observer,
                description="The Manning friction number applying to the current reach."
            ),
            base.Input(
                "BankSlope",
                (attrib.Class(list[float]), attrib.Unit("1"), attrib.Scales("space/base_geometry")),
                self.default_observer,
                description="The slope of the reach."
            ),
            base.Input(
                "Width",
                (attrib.Class(list[float]), attrib.Unit("m"), attrib.Scales("space/base_geometry")),
                self.default_observer,
                description="The width of the reach (undocumented by the module)."
            ),
            base.Input(
                "Shape",
                (
                    attrib.Class(list[str]),
                    attrib.Unit(None),
                    attrib.Scales("space/base_geometry"),
                    attrib.InList(("TriangularReach", "RectangularReach", "SWATReachType"))
                ),
                self.default_observer,
                description="The shape of the current reach."
            ),
            base.Input(
                "BulkDensity",
                (attrib.Class(list[float]), attrib.Unit("kg/m³"), attrib.Scales("space/base_geometry")),
                self.default_observer,
                description="The mass density of the reach sediment."
            ),
            base.Input(
                "Porosity",
                (attrib.Class(list[float]), attrib.Unit("m³/m³"), attrib.Scales("space/base_geometry")),
                self.default_observer,
                description="The porosity of the reach sediment."
            ),
            base.Input(
                "OrganicContent",
                (attrib.Class(list[float]), attrib.Unit("g/g"), attrib.Scales("space/base_geometry")),
                self.default_observer,
                description="The amount of organic material in the sediment of the reach."
            ),
            base.Input(
                "SedimentDepth1stLayer",
                (attrib.Class(list[float]), attrib.Unit("m"), attrib.Scales("space/base_geometry")),
                self.default_observer,
                description="The depth of the first layer of sediment."
            ),
            base.Input(
                "SedimentDepth2ndLayer",
                (attrib.Class(list[float]), attrib.Unit("m"), attrib.Scales("space/base_geometry")),
                self.default_observer,
                description="The depth of the second layer of sediment."
            )
        ])
        self._outputs = base.OutputContainer(self, [
            base.Output(
                "PEC_SW",
                store,
                self,
                {"data_type": np.float, "scales": "time/hour, space/base_geometry"},
                "The modelled concentration in the water phase.",
                {
                    "type": np.ndarray,
                    "shape": ("the number of simulated hours", "the number of simulated reaches"),
                    "chunks": "for fast retrieval of time series",
                    "unit": "mg/m³"
                }
            ),
            base.Output(
                "MASS_SW",
                store,
                self,
                {"data_type": np.float, "scales": "time/hour, space/base_geometry"},
                "The modelled substance mass in the water phase.",
                {
                    "type": np.ndarray,
                    "shape": ("the number of simulated hours", "the number of simulated reaches"),
                    "chunks": "for fast retrieval of time series",
                    "unit": "mg"
                }
            ),
            base.Output(
                "MASS_SED",
                store,
                self,
                {"data_type": np.float, "scales": "time/hour, space/base_geometry"},
                "The modelled substance mass in sediment.",
                {
                    "type": np.ndarray,
                    "shape": ("the number of simulated hours", "the number of simulated reaches"),
                    "chunks": "for fast retrieval of time series",
                    "unit": "mg"
                }
            ),
            base.Output(
                "MASS_SED_DEEP",
                store,
                self,
                {"data_type": np.float, "scales": "time/hour, space/base_geometry"},
                "The modelled substance mass in deep sediment.",
                {
                    "type": np.ndarray,
                    "shape": ("the number of simulated hours", "the number of simulated reaches"),
                    "chunks": "for fast retrieval of time series",
                    "unit": "mg"
                }
            ),
            base.Output(
                "PEC_SED",
                store,
                self,
                {"data_type": np.float, "scales": "time/hour, space/base_geometry"},
                "The modelled concentration in sediment.",
                {
                    "type": np.ndarray,
                    "shape": ("the number of simulated hours", "the number of simulated reaches"),
                    "chunks": "for fast retrieval of time series",
                    "unit": "mg/m³"
                }
            ),
            base.Output(
                "Reaches",
                store,
                self,
                {"scales": "space/reach"},
                "The numerical identifiers of the reaches in the order of the other outputs.",
                {"type": "list[int]"}
            )
        ])
        self._begin = None
        self._timeString = None
        return

    def run(self):
        """
        Runs the component.
        :return: Nothing.
        """
        self.default_observer.write_message(2, "Component relies on insensible high precision of z-coordinate")
        project_name = "h1"
        processing_path = self.inputs["ProcessingPath"].read().values
        project_path = os.path.join(processing_path, project_name)
        os.makedirs(project_path)
        self.prepare_reaches_and_drift_deposition(os.path.join(project_path, "HydroList.csv"),
                                                  os.path.join(project_path, "ReachList.csv"),
                                                  os.path.join(project_path, "SprayDriftList.csv"),
                                                  os.path.join(project_path, "DrainageList.csv"))
        self.prepare_catchment_list(os.path.join(project_path, "CatchmentList.csv"))
        self.prepare_project_list(processing_path, project_name)
        self.prepare_substance_list(os.path.join(project_path, "SubstanceList.csv"))
        self.run_project(processing_path, project_name)
        self.read_outputs(os.path.join(project_path, "h1_reaches.h5"))
        return

    def prepare_project_list(self, processing_path, project_name):
        """
        Prepares th project list.
        :param processing_path: The working directory of the module.
        :param project_name: The name of the module project.
        :return: Nothing.
        """
        project_list_file = os.path.join(processing_path, project_name + ".csv")
        with open(project_list_file, "w") as f:
            # noinspection SpellCheckingInspection
            f.write(
                "key,fpath,database,begin,end,t_input_start,timestep_input,timestep_simulation,timestep_output," +
                "aggregation_output,withHydro,write_summary,substance,simulation,preprocessing,postprocessing," +
                "MASS_SW,MASS_SED,MASS_SED_DEEP,PEC_SW,PEC_SED,threshold_sw,threshold_sed\n")
            f.write(project_name + ",")  # key
            f.write(processing_path + ",")  # file path
            f.write("csv,")  # database
            f.write(self._begin.strftime("%Y-%m-%dT%H:%M") + ",")  # begin
            f.write(self._timeString + ",")  # end
            f.write(self._begin.strftime("%Y-%m-%dT%H:%M") + ",")  # t_input_start
            f.write("1H,")  # time step input
            f.write("1Min,")  # time step simulation
            f.write("1H,")  # time step output
            f.write("max,")  # aggregation_output
            f.write("FALSE,")  # withHydro
            f.write("FALSE,")  # writeSummary
            f.write("CMP_A,")  # substance
            f.write("TRUE,")  # simulation
            f.write("FALSE,")  # pre-processing
            f.write("FALSE,")  # postprocessing
            f.write("TRUE,")  # *MASS_SW
            f.write("TRUE,")  # *MASS_SED
            f.write("TRUE,")  # *MASS_SED_DEEP
            f.write("TRUE,")  # PEC_SW
            f.write("TRUE,")  # PEC_SED
            f.write(str(self.inputs["ThresholdSW"].read().values) + ",")  # threshold_sw
            f.write(str(self.inputs["ThresholdSediment"].read().values) + "\n")  # threshold_sed
        return

    def prepare_reaches_and_drift_deposition(self, reaches_file, reach_list_file, spray_drift_file,drainage_list_file):
        """
        Prepares the reaches and drift deposition inputs.
        :param reaches_file: The file path of the reach file.
        :param reach_list_file: The file path of the reach list file.
        :param spray_drift_file: The file path of the spray-drift file.
        :return: Nothing.
        """
        hydrography = self.inputs["Hydrography"].read().values
        reaches_hydrology = self.inputs["ReachesHydrology"].read().values
        self._begin = self.inputs["TimeSeriesStart"].read().values
        number_time_steps = self.inputs["WaterDischarge"].describe()["shape"][0]
        reaches_drift = self.inputs["DriftDeposition"].describe()["element_names"][1].get_values()
        reaches_sorted = [int(r[1:]) for r in sorted([f"r{r}" for r in reaches_hydrology])]
        hydrography_geometries = self.inputs["HydrographyGeometries"].read()
        hydrography_reaches = hydrography_geometries.element_names[0].get_values()
        downstream_reaches = self.inputs["DownstreamReach"].read().values
        initial_depths = self.inputs["InitialDepth"].read().values
        manning = self.inputs["Manning"].read().values
        bank_slopes = self.inputs["BankSlope"].read().values
        widths = self.inputs["Width"].read().values
        shapes = self.inputs["Shape"].read().values
        bulk_densities = self.inputs["BulkDensity"].read().values
        porosity = self.inputs["Porosity"].read().values
        organic_contents = self.inputs["OrganicContent"].read().values
        depths_sediment_1 = self.inputs["SedimentDepth1stLayer"].read().values
        depths_sediment_2 = self.inputs["SedimentDepth2ndLayer"].read().values
        self.outputs["Reaches"].set_values(reaches_sorted, element_names=(self.outputs["Reaches"],))
        #reaches_drift = self.inputs["ReachesDrift"].read().values
        driver = ogr.GetDriverByName("ESRI Shapefile")
        data_source = driver.Open(hydrography, 0)
        layer = data_source.GetLayer()
        reaches_sorted = [int(r[1:]) for r in sorted(["r" + str(r) for r in reaches_hydrology])]
        #self.outputs["Reaches"].set_values(reaches_sorted)
        with open(reach_list_file, "w") as f:
            # noinspection SpellCheckingInspection
            f.write(
                "key,x,y,z,downstream,initial_depth,manning_n,bankslope,bottomwidth,floodplainslope,shape,dens," +
                "porosity,oc,depth_sed,depth_sed_deep\n")
            with open(reaches_file, "w") as f2:
                f2.write("key,time,volume,flow,area\n")
                with open(spray_drift_file, "w") as f3:
                    f3.write("key,substance,time,rate\n")
                    for reach in reaches_sorted:
                        hydrography_index = hydrography_reaches.index(reach)
                        geom = ogr.CreateGeometryFromWkb(hydrography_geometries.values[hydrography_index])
                        coord = geom.GetPoint(0)
                        downstream = downstream_reaches[hydrography_index]
                        f.write(f"r{reach},")
                        f.write(f"{round(coord[0], 2)},")
                        f.write(f"{round(coord[1], 2)},")
                        f.write(f"{round(coord[2], 8)},")
                        f.write(f"{'' if downstream == 'Outlet' else 'r'}{downstream},")
                        f.write(f"{initial_depths[hydrography_index]},")
                        f.write(f"{manning[hydrography_index]},")
                        f.write(f"{bank_slopes[hydrography_index]},")
                        f.write(f"{widths[hydrography_index]},")
                        f.write("200,")  # floodplain
                        f.write(f"{shapes[hydrography_index]},")
                        f.write(f"{bulk_densities[hydrography_index]},")
                        f.write(f"{porosity[hydrography_index]},")
                        f.write(f"{organic_contents[hydrography_index]},")
                        f.write(f"{depths_sediment_1[hydrography_index]},")
                        f.write(f"{depths_sediment_2[hydrography_index]}\n")
                        i = int(np.where(reaches_hydrology == reach)[0])
                        discharge = self.inputs["WaterDischarge"].read(slices=(slice(number_time_steps), i)).values
                        volume = self.inputs["WaterVolume"].read(slices=(slice(number_time_steps), i)).values
                        area = self.inputs["WetSurfaceArea"].read(slices=(slice(number_time_steps), i)).values
                        j = int(np.where(reaches_drift == reach)[0])
                        drift_deposition = self.inputs["DriftDeposition"].read(
                            slices=(slice(int(number_time_steps / 24)), j)).values
                        for t in range(number_time_steps):
                            self._timeString = (self._begin + datetime.timedelta(hours=t)).strftime(
                                "%Y-%m-%dT%H:%M")
                            f2.write(f"r{reach},")
                            f2.write(f"{self._timeString},")
                            f2.write(f"{round(float(volume[t]), 2)},")
                            f2.write(f"{round(float(discharge[t]), 2)},")
                            f2.write(f"{round(float(area[t]), 2)}\n")
                            if t % 24 == 11:
                                drift_deposition_value = drift_deposition[int((t - 11) / 24)]
                                if drift_deposition_value > 0:
                                    f3.write(f"r{reach},")
                                    f3.write("CMP_A,")
                                    f3.write(f"{self._timeString},")
                                    f3.write(f"{format(float(drift_deposition_value), 'f')}")
                                    f3.write("\n")
                    with open(drainage_list_file, "w") as f4:
                        f4.write("key,substance,time,rate\n")
                        for reach in reaches_sorted:
                            layer.SetAttributeFilter("key = '{}'".format(reach))
                            for feature in layer:
                                geom = feature.GetGeometryRef()
                                coord = geom.GetPoint(0)
                                downstream = feature.GetField("downstream")
                                f.write("r" + str(reach) + ",")
                                f.write(str(round(coord[0], 2)) + ",")  # x
                                f.write(str(round(coord[1], 2)) + ",")  # y
                                f.write(str(round(coord[2], 8)) + ",")  # z
                                f.write(("" if downstream == "Outlet" else "r") + downstream + ',')
                                f.write(str(feature.GetField("initial_de")) + ",")
                                f.write(str(feature.GetField("manning_n")) + ",")
                                # noinspection SpellCheckingInspection
                                f.write(str(feature.GetField("bankslope")) + ",")
                                f.write(str(feature.GetField("width")) + ",")
                                f.write("200,")  # floodplain
                                f.write(feature.GetField("shape_1") + ",")
                                f.write(str(feature.GetField("dens")) + ",")
                                f.write(str(feature.GetField("porosity")) + ",")
                                f.write(str(feature.GetField("oc")) + ",")
                                f.write(str(feature.GetField("depth_sed")) + ",")
                                f.write(str(feature.GetField("depth_sed_")) + "\n")
                                i = int(np.where(reaches_hydrology == reach)[0])
                                discharge = self.inputs["WaterDischarge"].read(slices=(slice(number_time_steps), i)).values
                                volume = self.inputs["WaterVolume"].read(slices=(slice(number_time_steps), i)).values
                                area = self.inputs["WetSurfaceArea"].read(slices=(slice(number_time_steps), i)).values
                                j = int(np.where(reaches_drift == reach)[0])
                                drift_deposition = self.inputs["DriftDeposition"].read(
                                    slices=(slice(int(number_time_steps / 24)), j)).values
                                for t in range(number_time_steps):
                                    self._timeString = (self._begin + datetime.timedelta(hours=t)).strftime(
                                        "%Y-%m-%dT%H:%M")
                                    f2.write("r" + str(reach) + ",")
                                    f2.write(self._timeString + ",")
                                    f2.write(str(round(float(volume[t]), 2)) + ",")
                                    f2.write(str(round(float(discharge[t]), 2)) + ",")
                                    f2.write(str(round(float(area[t]), 2)) + "\n")
                                    if t % 24 == 11:
                                        drift_deposition_value = drift_deposition[int((t - 11) / 24)]
                                        if drift_deposition_value > 0:
                                            f3.write("r" + str(reach) + ",")
                                            f3.write("CMP_A,")
                                            f3.write(self._timeString + ",")
                                            f3.write("{:f}".format(float(drift_deposition_value)))
                                            f3.write("\n")
                                            
                                                                # get drift results
                                DrainageLoad = self.inputs["DrainageLoad"].read( slices=(slice(number_time_steps), i)).values
                    
                                # write drainage list
                                for t in range(number_time_steps):
                                    self._timeString = (self._begin + datetime.timedelta(hours=t)).strftime(
                                        "%Y-%m-%dT%H:%M")
                                    if float(DrainageLoad[t]) > 0:
                                        f4.write("r" + str(reach) + ",")
                                        f4.write("CMP_A,")
                                        f4.write(self._timeString + ",")
                                        f4.write("%.2f"%float(DrainageLoad[t]))
                                        f4.write("\n")

                                        
                        layer.ResetReading()
        return

    def prepare_catchment_list(self, catchment_file):
        """
        Prepares the catchment list.
        :param catchment_file: The file path of the catchment list.
        :return: Nothing.
        """
        shutil.copyfile(self.inputs["Catchment"].read().values, catchment_file)
        return

    def run_project(self, processing_path, project_name):
        """
        Runs the module project.
        :param processing_path: The working directory for the module.
        :param project_name: The name of the module project.
        :return: Nothing.
        """
        python_exe = os.path.join(os.path.dirname(__file__), "module", "bin", "python", "python.exe")
        python_script = os.path.join(os.path.dirname(__file__), "module", "bin", "main.py")
        # noinspection SpellCheckingInspection
        base.run_process(
            (python_exe, python_script, "--folder", processing_path, "--runlist", project_name),
            None,
            self.default_observer,
            {"HOMEPATH": processing_path}
        )
        return

    def prepare_substance_list(self, substance_list_file):
        """
        Prepares the substance list.
        :param substance_list_file: The file path of the substance list.
        :return: Nothing.
        """
        with open(substance_list_file, "w") as f:
            # noinspection SpellCheckingInspection
            f.write("key,molarmass,DT50sw,DT50sed,KOC,Temp0,Q10,plantuptake,QFAC\n")
            # noinspection SpellCheckingInspection
            f.write("CMP_A,{},{},{},{},{},{},{},{}\n".format(
                self.inputs["MolarMass"].read().values,
                self.inputs["DT50sw"].read().values,
                self.inputs["DT50sed"].read().values,
                self.inputs["KOC"].read().values,
                self.inputs["Temp0"].read().values,
                self.inputs["Q10"].read().values,
                self.inputs["PlantUptake"].read().values,
                self.inputs["QFac"].read().values
            ))
        return

    def read_outputs(self, output_file):
        """
        Reads the module outputs into the landscape model.
        :param output_file: The file path of the module output file.
        :return: Nothing.
        """
        with h5py.File(output_file) as f:
            for variable in [
                ("PEC_SW", "mg/m³"),
                ("MASS_SW", "mg"),
                ("MASS_SED", "mg"),
                ("MASS_SED_DEEP", "mg"),
                ("PEC_SED", "mg/m³")
            ]:
                data = f["/" + variable[0]]
                self.outputs[variable[0]].set_values(
                    np.ndarray,
                    shape=data.shape,
                    chunks=(min(262144, data.shape[0]), 1),
                    unit=variable[1],
                    element_names=(None, self.outputs["Reaches"]),
                    offset=(self._begin, None)
                )
                for chunk in base.chunk_slices(data.shape, (min(262144, data.shape[0]), 1)):
                    self.outputs[variable[0]].set_values(data[chunk], slices=chunk, create=False, calculate_max=True)
        return
