"""Component for the Steps environmental fate module."""
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
    The component encapsulating the StepsRiverNetwork environmental fate module. StepsRiverNetwork simulates in-stream
    environmental fate processes of pesticides for an entire river network of a catchment. Environmental fate processes
    are calculated for each single reach, and transport across the entire catchment is reported in an explicit timestep
    of one hour.
    """
    # RELEASES
    VERSION = base.VersionCollection(
        base.VersionInfo("2.1.7", "2023-09-18"),
        base.VersionInfo("2.1.6", "2023-09-13"),
        base.VersionInfo("2.1.5", "2023-09-12"),
        base.VersionInfo("2.1.4", "2023-09-11"),
        base.VersionInfo("2.1.3", "2022-03-03"),
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
    VERSION.changed("2.1.3", "Mitigated weak code warnings")
    VERSION.added("2.1.4", "Information on runtime environment")
    VERSION.changed("2.1.5", "Extended module information for Python runtime environment")
    VERSION.added("2.1.5", "Creation of repository info during documentation")
    VERSION.added("2.1.5", "Repository info, changelog, contributing note and readme to module")
    VERSION.added("2.1.5", "Repository info to Python runtime environment")
    VERSION.added("2.1.6", "Scales attribute to global inputs")
    VERSION.fixed("2.1.6", "Spatial scale of outputs")
    VERSION.changed("2.1.7", "Updated component description")
    VERSION.changed("2.1.7", "Updated input descriptions and removed stub descriptions")
    VERSION.added("2.1.7", "Runtime warnings and notes regarding status of component and documentation")

    def __init__(self, name, observer, store):
        """
        Initializes a StepsRiverNetwork component.

        Args:
            name: The name of the component.
            observer: The default observer of the component.
            store: The default store of the component.
        """
        super(StepsRiverNetwork, self).__init__(name, observer, store)
        self._module = base.Module(
            "River network version of STEPS1234",
            "0.93",
            "module",
            r"module\documentation\html\index.html",
            base.Module(
                "Python",
                "3.7.4",
                "module/bin/python",
                "module/bin/python/Doc/python374.chm",
                None,
                True,
                "module/bin/python/NEWS.txt"
            )
        )
        # noinspection SpellCheckingInspection
        self._inputs = base.InputContainer(self, [
            base.Input(
                "ProcessingPath",
                (attrib.Class(str), attrib.Unit(None), attrib.Scales("global")),
                self.default_observer,
                description="The working directory for the module. It is used for all files prepared as module inputs "
                            "or generated as (temporary) module outputs."
            ),
            base.Input(
                "Catchment",
                (attrib.Class(str), attrib.Unit(None), attrib.Scales("global")),
                self.default_observer,
                description="A file path to a CSV file detailing the hydrographic properties of the entire catchment "
                            "depicted by a hydrographic scenario. This file is usually provided by the scenario "
                            "developer (if usage of StepsRiverNetwork is supported by the scenario) and is made "
                            "available as a project macro. See the module documentation for details on the format."
            ),
            base.Input(
                "WaterDischarge",
                (attrib.Class(np.ndarray), attrib.Unit("m³/d"), attrib.Scales("time/hour, space/reach")),
                self.default_observer
            ),
            base.Input(
                "TimeSeriesStart",
                (attrib.Class(datetime.datetime), attrib.Unit(None), attrib.Scales("global")),
                self.default_observer,
                description="The first time step for which input data is provided. This is also the time step of where "
                            "the StepsRiverNetwork simulation starts. This input will be removed in a future version "
                            "of the `StepsRiverNetwork` component."
            ),
            base.Input(
                "ReachesHydrology",
                (attrib.Class(np.ndarray, 1), attrib.Unit(None, 1), attrib.Scales("space/reach", 1)),
                self.default_observer,
                description="The numeric identifiers for individual reaches (in the order of the hydrological inputs). "
                            "his input will be removed in a future version of the `StepsRiverNetwork` component."
            ),
            base.Input(
                "WaterVolume",
                (attrib.Class(np.ndarray), attrib.Unit("m³"), attrib.Scales("time/hour, space/reach")),
                self.default_observer
            ),
            base.Input(
                "WetSurfaceArea",
                (attrib.Class(np.ndarray), attrib.Unit("m²"), attrib.Scales("time/hour, space/reach")),
                self.default_observer
            ),
            base.Input(
                "DriftDeposition",
                (attrib.Class(np.ndarray), attrib.Unit("mg/m²"), attrib.Scales("time/day, space/reach")),
                self.default_observer
            ),
            base.Input(
                "MolarMass",
                (attrib.Class(float), attrib.Unit("g/mol"), attrib.Scales("global")),
                self.default_observer
            ),
            base.Input(
                "DT50sw", (attrib.Class(float), attrib.Unit("d"), attrib.Scales("global")), self.default_observer),
            base.Input(
                "DT50sed", (attrib.Class(float), attrib.Unit("d"), attrib.Scales("global")), self.default_observer),
            base.Input(
                "KOC", (attrib.Class(float), attrib.Unit("l/kg"), attrib.Scales("global")), self.default_observer),
            base.Input(
                "Temp0", (attrib.Class(float), attrib.Unit("°C"), attrib.Scales("global")), self.default_observer),
            base.Input(
                "Q10", (attrib.Class(float, 1), attrib.Unit("1", 1), attrib.Scales("global")), self.default_observer),
            base.Input(
                "PlantUptake", (attrib.Class(float), attrib.Unit("1"), attrib.Scales("global")), self.default_observer),
            base.Input(
                "QFac", (attrib.Class(float), attrib.Unit("1"), attrib.Scales("global")), self.default_observer),
            base.Input(
                "ThresholdSW",
                (attrib.Class(float), attrib.Unit("mg/m³"), attrib.Scales("global")),
                self.default_observer
            ),
            base.Input(
                "ThresholdSediment",
                (attrib.Class(float), attrib.Unit("mg/kg"), attrib.Scales("global")),
                self.default_observer
            ),
            base.Input(
                "HydrographyGeometries",
                (attrib.Class(list[bytes]), attrib.Unit(None), attrib.Scales("space/reach")),
                self.default_observer,
                description="The geometries of individual water body segments (reaches) in WKB representation. This "
                            "input will be removed in a future version of the `StepsRiverNetwork` component."
            ),
            base.Input(
                "DownstreamReach",
                (attrib.Class(list[str]), attrib.Unit(None), attrib.Scales("space/reach")),
                self.default_observer
            ),
            base.Input(
                "InitialDepth",
                (attrib.Class(list[float]), attrib.Unit("m"), attrib.Scales("space/reach")),
                self.default_observer
            ),
            base.Input(
                "Manning",
                (attrib.Class(list[float]), attrib.Unit("1"), attrib.Scales("space/reach")),
                self.default_observer
            ),
            base.Input(
                "BankSlope",
                (attrib.Class(list[float]), attrib.Unit("1"), attrib.Scales("space/reach")),
                self.default_observer
            ),
            base.Input(
                "Width",
                (attrib.Class(list[float]), attrib.Unit("m"), attrib.Scales("space/reach")),
                self.default_observer
            ),
            base.Input(
                "Shape",
                (
                    attrib.Class(list[str]),
                    attrib.Unit(None),
                    attrib.Scales("space/reach"),
                    attrib.InList(("TriangularReach", "RectangularReach", "SWATReachType"))
                ),
                self.default_observer
            ),
            base.Input(
                "BulkDensity",
                (attrib.Class(list[float]), attrib.Unit("kg/m³"), attrib.Scales("space/reach")),
                self.default_observer
            ),
            base.Input(
                "Porosity",
                (attrib.Class(list[float]), attrib.Unit("m³/m³"), attrib.Scales("space/reach")),
                self.default_observer
            ),
            base.Input(
                "OrganicContent",
                (attrib.Class(list[float]), attrib.Unit("g/g"), attrib.Scales("space/reach")),
                self.default_observer
            ),
            base.Input(
                "SedimentDepth1stLayer",
                (attrib.Class(list[float]), attrib.Unit("m"), attrib.Scales("space/reach")),
                self.default_observer
            ),
            base.Input(
                "SedimentDepth2ndLayer",
                (attrib.Class(list[float]), attrib.Unit("m"), attrib.Scales("space/reach")),
                self.default_observer
            )
        ])
        self._outputs = base.OutputContainer(self, [
            base.Output(
                "PEC_SW",
                store,
                self,
                {"data_type": np.float, "scales": "time/hour, space/reach"},
                "The modelled concentration in the water phase.",
                {
                    "type": np.ndarray,
                    "shape": ("the number of simulated hours", "the number of simulated reaches"),
                    "chunks": "for fast retrieval of time series",
                    "unit": "mg/m³",
                    "element_names": (None, "as specified by the `Reaches` output"),
                    "offset": ("the value of the `TimeSeriesStart` input", None),
                    "geometries": (None, "as specified by the `ReachesGeometries` output")
                }
            ),
            base.Output(
                "MASS_SW",
                store,
                self,
                {"data_type": np.float, "scales": "time/hour, space/reach"},
                "The modelled substance mass in the water phase.",
                {
                    "type": np.ndarray,
                    "shape": ("the number of simulated hours", "the number of simulated reaches"),
                    "chunks": "for fast retrieval of time series",
                    "unit": "mg",
                    "element_names": (None, "as specified by the `Reaches` output"),
                    "offset": ("the value of the `TimeSeriesStart` input", None),
                    "geometries": (None, "as specified by the `ReachesGeometries` output")
                }
            ),
            base.Output(
                "MASS_SED",
                store,
                self,
                {"data_type": np.float, "scales": "time/hour, space/reach"},
                "The modelled substance mass in sediment.",
                {
                    "type": np.ndarray,
                    "shape": ("the number of simulated hours", "the number of simulated reaches"),
                    "chunks": "for fast retrieval of time series",
                    "unit": "mg",
                    "element_names": (None, "as specified by the `Reaches` output"),
                    "offset": ("the value of the `TimeSeriesStart` input", None),
                    "geometries": (None, "as specified by the `ReachesGeometries` output")
                }
            ),
            base.Output(
                "MASS_SED_DEEP",
                store,
                self,
                {"data_type": np.float, "scales": "time/hour, space/reach"},
                "The modelled substance mass in deep sediment.",
                {
                    "type": np.ndarray,
                    "shape": ("the number of simulated hours", "the number of simulated reaches"),
                    "chunks": "for fast retrieval of time series",
                    "unit": "mg",
                    "element_names": (None, "as specified by the `Reaches` output"),
                    "offset": ("the value of the `TimeSeriesStart` input", None),
                    "geometries": (None, "as specified by the `ReachesGeometries` output")
                }
            ),
            base.Output(
                "PEC_SED",
                store,
                self,
                {"data_type": np.float, "scales": "time/hour, space/reach"},
                "The modelled concentration in sediment.",
                {
                    "type": np.ndarray,
                    "shape": ("the number of simulated hours", "the number of simulated reaches"),
                    "chunks": "for fast retrieval of time series",
                    "unit": "mg/m³",
                    "element_names": (None, "as specified by the `Reaches` output"),
                    "offset": ("the value of the `TimeSeriesStart` input", None),
                    "geometries": (None, "as specified by the `ReachesGeometries` output")
                }
            ),
            base.Output(
                "Reaches",
                store,
                self,
                {"scales": "space/reach"},
                "The numerical identifiers of the reaches in the order of the other outputs.",
                {
                    "type": list[int],
                    "element_names": ("as specified by the output itself",),
                    "geometries": ("as specified by the `ReachesGeometries` output",)
                }
            ),
            base.Output(
                "ReachesGeometries",
                store,
                self,
                {"scales": "space/reach"},
                "The geometry of the reaches in Well-Known-Byte notation in the order of the reaches as specified by "
                "the `Reaches` output.",
                {
                    "type": list[bytes],
                    "element_names": ("as specified by the `Reaches` output",),
                    "geometries": ("as specified by output itself",)
                }
            )
        ])
        self._begin = None
        self._timeString = None
        if self.default_observer:
            self.default_observer.write_message(
                2,
                "StepsRiverNetwork currently does not check the identity of reaches",
                "Make sure that inputs of scale space/reach retrieve data in the same reach-order"
            )
            self.default_observer.write_message(
                3,
                "The TimeSeriesStart input will be removed in a future version of the StepsRiverNetwork component",
                "The time offset will be retrieved from the metadata of the WaterDischarge input"
            )
            self.default_observer.write_message(
                3,
                "The ReachesHydrology input will be removed in a future version of the StepsRiverNetwork component",
                "The reach names will be retrieved from the metadata of the WaterDischarge input"
            )
            self.default_observer.write_message(
                3,
                "The HydrographyGeometries input will be removed in a future version of the StepsRiverNetwork "
                "component",
                "The reach geometries will be retrieved from the metadata of the WaterDischarge input"
            )

    def run(self):
        """
        Runs the component.

        Returns:
            Nothing.
        """
        self.default_observer.write_message(2, "Component relies on insensible high precision of z-coordinate")
        project_name = "h1"
        processing_path = self.inputs["ProcessingPath"].read().values
        project_path = os.path.join(processing_path, project_name)
        os.makedirs(project_path)
        self.prepare_reaches_and_drift_deposition(
            os.path.join(project_path, "HydroList.csv"),
            os.path.join(project_path, "ReachList.csv"),
            os.path.join(project_path, "SprayDriftList.csv")
        )
        self.prepare_catchment_list(os.path.join(project_path, "CatchmentList.csv"))
        self.prepare_project_list(processing_path, project_name)
        self.prepare_substance_list(os.path.join(project_path, "SubstanceList.csv"))
        self.run_project(processing_path, project_name)
        self.read_outputs(os.path.join(project_path, "h1_reaches.h5"))

    def prepare_project_list(self, processing_path, project_name):
        """
        Prepares th project list.

        Args:
            processing_path: The working directory of the module.
            project_name: The name of the module project.

        Returns:
            Nothing.
        """
        project_list_file = os.path.join(processing_path, f"{project_name}.csv")
        with open(project_list_file, "w") as f:
            # noinspection SpellCheckingInspection
            f.write(
                "key,fpath,database,begin,end,t_input_start,timestep_input,timestep_simulation,timestep_output,"
                "aggregation_output,withHydro,write_summary,substance,simulation,preprocessing,postprocessing,"
                "MASS_SW,MASS_SED,MASS_SED_DEEP,PEC_SW,PEC_SED,threshold_sw,threshold_sed\n"
            )
            f.write(f"{project_name},")
            f.write(f"{processing_path},")
            f.write("csv,")  # database
            f.write(f"{self._begin.strftime('%Y-%m-%dT%H:%M')},")
            f.write(f"{self._timeString},")
            f.write(f"{self._begin.strftime('%Y-%m-%dT%H:%M')},")
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
            f.write(f"{self.inputs['ThresholdSW'].read().values},")
            f.write(f"{self.inputs['ThresholdSediment'].read().values}\n")

    def prepare_reaches_and_drift_deposition(self, reaches_file, reach_list_file, spray_drift_file):
        """
        Prepares the reaches and drift deposition inputs.

        Args:
            reaches_file: The file path of the reach file.
            reach_list_file: The file path of the reach list file.
            spray_drift_file: The file path of the spray-drift file.

        Returns:
            Nothing.
        """
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
        geometries_sorted = [hydrography_geometries.values[hydrography_reaches.index(x)] for x in reaches_sorted]
        self.outputs["Reaches"].set_values(
            reaches_sorted, element_names=(self.outputs["Reaches"],), geometries=(self.outputs["ReachesGeometries"],))
        self.outputs["ReachesGeometries"].set_values(
            geometries_sorted,
            element_names=(self.outputs["Reaches"],),
            geometries=(self.outputs["ReachesGeometries"],)
        )
        with open(reach_list_file, "w") as f:
            # noinspection SpellCheckingInspection
            f.write(
                "key,x,y,z,downstream,initial_depth,manning_n,bankslope,bottomwidth,floodplainslope,shape,dens,"
                "porosity,oc,depth_sed,depth_sed_deep\n"
            )
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
                        f.write(f"{round(coord[0], 2)},{round(coord[1], 2)},{round(coord[2], 8)},")
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

    def prepare_catchment_list(self, catchment_file):
        """
        Prepares the catchment list.

        Args:
            catchment_file: The file path of the catchment list.

        Returns:
            Nothing.
        """
        shutil.copyfile(self.inputs["Catchment"].read().values, catchment_file)

    def run_project(self, processing_path, project_name):
        """
        Runs the module project.

        Args:
            processing_path: The working directory for the module.
            project_name: The name of the module project.

        Returns:
            Nothing.
        """
        python_exe = os.path.join(os.path.dirname(__file__), "module", "bin", "python", "python.exe")
        python_script = os.path.join(os.path.dirname(__file__), "module", "bin", "main.py")
        # noinspection SpellCheckingInspection
        base.run_process(
            (python_exe, python_script, "--folder", processing_path, "--runlist", project_name),
            processing_path,
            self.default_observer,
            {"HOMEPATH": processing_path}
        )

    # noinspection DuplicatedCode
    def prepare_substance_list(self, substance_list_file):
        """
        Prepares the substance list.

        Args:
            substance_list_file: The file path of the substance list.

        Returns:
            Nothing.
        """
        with open(substance_list_file, "w") as f:
            # noinspection SpellCheckingInspection
            f.write("key,molarmass,DT50sw,DT50sed,KOC,Temp0,Q10,plantuptake,QFAC\n")
            # noinspection SpellCheckingInspection
            f.write(
                f"CMP_A,{self.inputs['MolarMass'].read().values},{self.inputs['DT50sw'].read().values},"
                f"{self.inputs['DT50sed'].read().values},{self.inputs['KOC'].read().values},"
                f"{self.inputs['Temp0'].read().values},{self.inputs['Q10'].read().values},"
                f"{self.inputs['PlantUptake'].read().values},{self.inputs['QFac'].read().values}\n"
            )

    def read_outputs(self, output_file):
        """
        Reads the module outputs into the landscape model.

        Args:
            output_file: The file path of the module output file.

        Returns:
            Nothing.
        """
        with h5py.File(output_file) as f:
            for variable in [
                ("PEC_SW", "mg/m³"),
                ("MASS_SW", "mg"),
                ("MASS_SED", "mg"),
                ("MASS_SED_DEEP", "mg"),
                ("PEC_SED", "mg/m³")
            ]:
                data = f[f"/{variable[0]}"]
                self.outputs[variable[0]].set_values(
                    np.ndarray,
                    shape=data.shape,
                    chunks=(min(262144, data.shape[0]), 1),
                    unit=variable[1],
                    element_names=(None, self.outputs["Reaches"]),
                    offset=(self._begin, None),
                    geometries=(None, self.outputs["ReachesGeometries"])
                )
                for chunk in base.chunk_slices(data.shape, (min(262144, data.shape[0]), 1)):
                    self.outputs[variable[0]].set_values(data[chunk], slices=chunk, create=False, calculate_max=True)
