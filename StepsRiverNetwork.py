"""
Component for the Steps efate module.
"""
import datetime
import h5py
import numpy as np
import ogr
import os
import shutil
import base
import attrib


class StepsRivernetwork(base.Component):
    """
    The component encapsulating the Steps efate module.
    """
    # RELEASES
    VERSION = base.VersionCollection(
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

    # CHANGELOG
    VERSION.added("1.2.3", "components.StepsRivernetwork component")
    VERSION.added("1.2.4", "components.CmfHydrology hydrology Python library only loaded when needed")
    VERSION.changed("1.2.36", "components.CmfHydrology replaced by components.StepsRiverNetwork")
    VERSION.changed("1.2.37", "components.StepsRiverNetwork module updated to version 0.9.2")
    VERSION.changed("1.2.37", "components.StepsRiverNetwork thresholds as inputs")
    VERSION.fixed("1.3.3", "increased numeric precision in components.StepsRivernetwork (preliminary)")
    VERSION.changed("1.3.7", "components.StepsRiverNetwork module updated to version 0.9.3")
    VERSION.fixed("1.3.10", "Spatial referencing in components.StepsRiverNetwork")
    VERSION.changed("1.3.16", "Substance parameterization in components.StepsRiverNetwork changed")
    VERSION.changed("1.3.24", "components.StepsRiverNetwork uses base function to call module")
    VERSION.fixed("1.3.26", "components.StepsRiverNetwork reach order")
    VERSION.changed("1.3.27", "components.StepsRiverNetwork specifies scales")
    VERSION.fixed("1.3.29", "Input slicing in components.StepsRiverNetwork")
    VERSION.changed("1.3.33", "components.StepsRiverNetwork checks input types strictly")
    VERSION.changed("1.3.33", "components.StepsRiverNetwork checks for physical units")
    VERSION.changed("1.3.33", "components.StepsRiverNetwork reports physical units to the data store")
    VERSION.changed("1.3.33", "components.StepsRiverNetwork checks for scales")
    VERSION.changed("1.3.35", "components.StepsRiverNetwork receives processing path as home path environment variable")
    VERSION.changed("2.0.0", "First independent release")
    VERSION.added("2.0.1", "Changelog and release history")

    def __init__(self, name, observer, store):
        super(StepsRivernetwork, self).__init__(name, observer, store)
        self._module = base.Module("River network version of STEPS1234", "0.93")
        # noinspection SpellCheckingInspection
        self._inputs = base.InputContainer(self, [
            base.Input(
                "ProcessingPath",
                (attrib.Class(str, 1), attrib.Unit(None, 1), attrib.Scales("global", 1)),
                self.default_observer
            ),
            base.Input(
                "Hydrography",
                (attrib.Class(str, 1), attrib.Unit(None, 1), attrib.Scales("global", 1)),
                self.default_observer
            ),
            base.Input(
                "Catchment",
                (attrib.Class(str, 1), attrib.Unit(None, 1), attrib.Scales("global", 1)),
                self.default_observer
            ),
            base.Input(
                "WaterDischarge",
                (
                    attrib.Class(np.ndarray, 1),
                    attrib.Unit("m³/d", 1),
                    attrib.Scales("time/hour, space/reach", 1)
                ),
                self.default_observer
            ),
            base.Input(
                "TimeseriesStart",
                (attrib.Class(datetime.datetime, 1), attrib.Unit(None, 1), attrib.Scales("global", 1)),
                self.default_observer
            ),
            base.Input(
                "ReachesHydrology",
                (attrib.Class(np.ndarray, 1), attrib.Unit(None, 1), attrib.Scales("space/reach", 1)),
                self.default_observer
            ),
            base.Input(
                "WaterVolume",
                (
                    attrib.Class(np.ndarray, 1),
                    attrib.Unit("m³", 1),
                    attrib.Scales("time/hour, space/reach", 1)
                ),
                self.default_observer
            ),
            base.Input(
                "WetSurfaceArea",
                (
                    attrib.Class(np.ndarray, 1),
                    attrib.Unit("m²", 1),
                    attrib.Scales("time/hour, space/reach", 1)
                ),
                self.default_observer
            ),
            base.Input(
                "DriftDeposition",
                (
                    attrib.Class(np.ndarray, 1),
                    attrib.Unit("mg/m²", 1),
                    attrib.Scales("time/day, space/reach", 1)
                ),
                self.default_observer
            ),
            base.Input(
                "ReachesDrift",
                (attrib.Class(np.ndarray, 1), attrib.Unit(None, 1), attrib.Scales("space/reach", 1)),
                self.default_observer
            ),
            base.Input("MolarMass", (attrib.Class(float, 1), attrib.Unit("g/mol", 1)), self.default_observer),
            base.Input("DT50sw", (attrib.Class(float, 1), attrib.Unit("d", 1)), self.default_observer),
            base.Input("DT50sed", (attrib.Class(float, 1), attrib.Unit("d", 1)), self.default_observer),
            base.Input("KOC", (attrib.Class(float, 1), attrib.Unit("l/kg", 1)), self.default_observer),
            base.Input("Temp0", (attrib.Class(float, 1), attrib.Unit("°C", 1)), self.default_observer),
            base.Input("Q10", (attrib.Class(float, 1), attrib.Unit("1", 1)), self.default_observer),
            base.Input("PlantUptake", (attrib.Class(float, 1), attrib.Unit("1", 1)), self.default_observer),
            base.Input("QFAC", (attrib.Class(float, 1), attrib.Unit("1", 1)), self.default_observer),
            base.Input(
                "ThresholdSW",
                (attrib.Class(float, 1), attrib.Unit("mg/m³", 1)),
                self.default_observer
            ),
            base.Input(
                "ThresholdSediment",
                (attrib.Class(float, 1), attrib.Unit("mg/kg", 1)),
                self.default_observer
            )
        ])
        self._outputs = base.OutputContainer(self, [
            base.Output("PEC_SW", store, self),
            base.Output("MASS_SW", store, self),
            base.Output("MASS_SED", store, self),
            base.Output("MASS_SED_DEEP", store, self),
            base.Output("PEC_SED", store, self),
            base.Output("Reaches", store, self)
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
                                                  os.path.join(project_path, "SprayDriftList.csv"))
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

    def prepare_reaches_and_drift_deposition(self, reaches_file, reach_list_file, spray_drift_file):
        """
        Prepares the reaches and drift deposition inputs.
        :param reaches_file: The file path of the reach file.
        :param reach_list_file: The file path of the reach list file.
        :param spray_drift_file: The file path of the spray-drift file.
        :return: Nothing.
        """
        hydrography = self.inputs["Hydrography"].read().values
        reaches_hydrology = self.inputs["ReachesHydrology"].read().values
        self._begin = self.inputs["TimeseriesStart"].read().values
        number_time_steps = self.inputs["WaterDischarge"].describe()["shape"][0]
        reaches_drift = self.inputs["ReachesDrift"].read().values
        driver = ogr.GetDriverByName("ESRI Shapefile")
        data_source = driver.Open(hydrography, 0)
        layer = data_source.GetLayer()
        reaches_sorted = [int(r[1:]) for r in sorted(["r" + str(r) for r in reaches_hydrology])]
        self.outputs["Reaches"].set_values(reaches_sorted, scales="space/reach")
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
                self.inputs["QFAC"].read().values
            ))
        return

    def read_outputs(self, output_file):
        """
        Reads the module outputs into the landscape model.
        :param output_file: The file path of the module output file.
        :return: Nothing.
        """
        with h5py.File(output_file, "r") as f:
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
                    data_type=np.float,
                    chunks=(min(262144, data.shape[0]), 1),
                    scales="time/hour, space/base_geometry",
                    unit=variable[1]
                )
                for chunk in base.chunk_slices(data.shape, (min(262144, data.shape[0]), 1)):
                    self.outputs[variable[0]].set_values(data[chunk], slices=chunk, create=False, calculate_max=True)
        return