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
    """The component encapsulating the Steps environmental fate module."""
    # RELEASES
    VERSION = base.VersionCollection(
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
                input basically represents the flow-lines used during preparation of the hydrology. The hydrography is
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
        self.prepare_reaches_and_drift_deposition(os.path.join(project_path, "HydroList.csv"),
                                                  os.path.join(project_path, "ReachList.csv"),
                                                  os.path.join(project_path, "SprayDriftList.csv"))
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
        hydrography = self.inputs["Hydrography"].read().values
        reaches_hydrology = self.inputs["ReachesHydrology"].read().values
        self._begin = self.inputs["TimeSeriesStart"].read().values
        number_time_steps = self.inputs["WaterDischarge"].describe()["shape"][0]
        reaches_drift = self.inputs["ReachesDrift"].read().values
        driver = ogr.GetDriverByName("ESRI Shapefile")
        data_source = driver.Open(hydrography, 0)
        layer = data_source.GetLayer()
        reaches_sorted = [int(r[1:]) for r in sorted([f"r{r}" for r in reaches_hydrology])]
        self.outputs["Reaches"].set_values(reaches_sorted)
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
                        layer.SetAttributeFilter(f"key = '{reach}'")
                        for feature in layer:
                            geom = feature.GetGeometryRef()
                            coord = geom.GetPoint(0)
                            downstream = feature.GetField("downstream")
                            f.write(f"r{reach},")
                            f.write(f"{round(coord[0], 2)},")
                            f.write(f"{round(coord[1], 2)},")
                            f.write(f"{round(coord[2], 8)},")
                            f.write(f"{'' if downstream == 'Outlet' else 'r'}{downstream},")
                            f.write(f"{feature.GetField('initial_de')},")
                            f.write(f"{feature.GetField('manning_n')},")
                            # noinspection SpellCheckingInspection
                            f.write(f"{feature.GetField('bankslope')},")
                            f.write(f"{feature.GetField('width')},")
                            f.write("200,")  # floodplain
                            f.write(f"{feature.GetField('shape_1')},")
                            f.write(f"{feature.GetField('dens')},")
                            f.write(f"{feature.GetField('porosity')},")
                            f.write(f"{feature.GetField('oc')},")
                            f.write(f"{feature.GetField('depth_sed')},")
                            f.write(f"{feature.GetField('depth_sed_')}\n")
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
                        layer.ResetReading()

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
            None,
            self.default_observer,
            {"HOMEPATH": processing_path}
        )

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
                    unit=variable[1]
                )
                for chunk in base.chunk_slices(data.shape, (min(262144, data.shape[0]), 1)):
                    self.outputs[variable[0]].set_values(data[chunk], slices=chunk, create=False, calculate_max=True)
