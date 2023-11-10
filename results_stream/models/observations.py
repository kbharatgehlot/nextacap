import json
import os
from pathlib import Path
import subprocess

from glob import glob
import numpy as np
import pandas as pd
import streamlit as st
import streamlit.components.v1 as components
from .gains import (
    Gains,
    get_all_stations_amplitude_stats,
    get_all_clusters_amplitude_stats,
    plot_gains_amplitude_stat_for_all_stations,
    plot_gains_amplitude_stat_for_all_clusters,
)
from .flags import Flags
from .utils import all_sky_map, convert_pdf_to_jpeg


# @st.cache_resource
def mock_jsons(size):
    js = [
        # "/home/users/chege/theleap/testdir/L254871_v4/calibrate_pb_spatial_deg2_t16_logs/params.json"
        # "/home/users/chege/theleap/testdir/L254871/logs/params.json",
        # "/home/users/chege/theleap/testdir/L254871/logs/params.json",
        # "/home/users/chege/theleap/testdir/L254871/logs/params.json",
        # "/home/users/chege/theleap/testdir/L254871_v2/logs/params.json",
        # "/home/users/chege/NCP/redshift1/L412388/logs/params2.json",
        "/home/users/chege/theleap/leap/test/params_node115.json"
    ]
    o = {}

    for i in range(size):
        n = "L2554" + str(i)
        o[n] = js[0]  # [np.random.randint(0, 1)]
    return o


all_calibration_stages = [
    "mpi_di",
    "bandpass",
    "mpi_dd",
]


class Processed:
    nights = mock_jsons(3)


class Config:
    def __init__(self, nextleap_config_file):
        self.nextleap_config_file = nextleap_config_file

        with open(self.nextleap_config_file, "r") as f:
            self.config = json.load(f)

        # tasks run
        self.tasks = self.config["tasks"]

        # observation name
        self.obsid = self.config["data"]["obsid"]

        # what calibration stages were done
        self.calibration_stages_done = list(
            set(self.tasks) & set(all_calibration_stages)
        )

        # nodes used
        self.working_nodes = self.config["cluster"]["nodes"].split(",")
        # master node
        self.master_node = self.config["cluster"]["masternode"]
        # data directory
        self.data_directory = (
            f"/net/{self.config['cluster']['masternode']}/{self.config['data']['path']}/results"
        )

        if "mpi_di" in self.tasks:
            self.di_solsdir = (
                f"{self.data_directory}/{self.config['mpi_di']['solsdir']}"
            )
            self.di_clusters_file = self.config["mpi_di"]["clusters_file"]
            self.di_solutions_file = f"{self.di_solsdir}/{self.obsid}.npy"
            _x = os.path.basename(self.di_clusters_file)
            self.di_effective_clusters = f"{self.di_solsdir}/eff_nr_{_x}.npy"

        if "mpi_dd" in self.tasks:
            self.dd_solsdir = (
                f"{self.data_directory}/{self.config['mpi_dd']['solsdir']}"
            )
            self.dd_clusters_file = self.config["mpi_dd"]["clusters_file"]
            self.dd_solutions_file = f"{self.dd_solsdir}/{self.obsid}.npy"
            _x = os.path.basename(self.dd_clusters_file)
            self.dd_effective_clusters = f"{self.dd_solsdir}/eff_nr_{_x}.npy"

        if "bandpass" in self.tasks:
            self.bandpass_solsdir = (
                f"{self.data_directory}/{self.config['bandpass']['solsdir']}"
            )
            self.bandpass_clusters_file = self.config["bandpass"]["clusters_file"]
            self.bandpass_solutions_file = f"{self.bandpass_solsdir}/{self.obsid}.npy"
            _x = os.path.basename(self.bandpass_clusters_file)
            self.bandpass_effective_clusters = (
                f"{self.bandpass_solsdir}/eff_nr_{_x}.npy"
            )

        if "ps" in self.tasks:
            self.pspipe_pdfs = sorted(
                [
                    path
                    for path in Path(
                        f"{self.data_directory}/{self.config['pspipe']['dir']}"
                    ).rglob("*.pdf")
                ]
            )

        if "wsclean" in self.tasks and "bandpass" in self.tasks:
            self.di_fits_image = glob(
                f"{self.data_directory}/{self.config['wsclean']['dir']}_DI_images/*image.fits"
            )
        else:
            self.di_fits_image = None

        if "wsclean" in self.tasks and "mpi_dd" in self.tasks:
            self.dd_fits_image = glob(
                f"{self.data_directory}/{self.config['wsclean']['dir']}_DD_images/*image.fits"
            )
        else:
            self.dd_fits_image = None

    def get_nextleap_report(self):
        return f"{self.config['logs_dir']}/report.html"

    def get_nextleap_trace(self):
        return f"{self.config['logs_dir']}/trace.html"

    def get_nextleap_timeline(self):
        return f"{self.config['logs_dir']}/timeline.html"


class ObsID(Config):
    def __init__(self, nextleap_config_file):
        super().__init__(nextleap_config_file)

    def __str__(self):
        return f"Obsid: {self.obsid} \nTasks: {self.tasks} \nData directory {self.data_directory} \nFull configuration {self.config})"

    def diagnostics_display(self):
        st.markdown(
            "<h2 style='text-align: center; color: velvet;'>NextLEAP 2023 processing cycle</h2>",
            unsafe_allow_html=True,
        )
        c1, c2 = st.columns([0.89, 0.1])
        c1.markdown(
            f"<h3 style='color: green;'>{self.obsid}</h3>",
            unsafe_allow_html=True,
        )
        if c2.button("Clear cache"):
            st.cache_data.clear()
        with st.container():
            with st.expander("Observation summary"):
                all_ms_files = self.all_ms_files(self.calibration_stages_done[-1])
                st.write(all_ms_files)
                # TODO: make this faster and formated
                # st.markdown(
                #     subprocess.check_output(
                #         f"msoverview in={all_ms_files[0]} verbose=T",
                #         stderr=subprocess.STDOUT,
                #         shell=True,
                #     ).decode("utf-8")
                # )

            with st.expander("NextLEAP configuration"):
                st.json(self.config)
            with st.expander("Nextleap report"):
                # html_file=self.get_nextleap_report()
                # self.display_html_file(html_file)
                trace_file = self.get_nextleap_trace()
                self.display_trace_file(trace_file)

                # timeline = self.get_nextleap_timeline()
                # self.display_html_file(timeline)
            with st.expander("RFI Occupancy"):
                for cal in self.calibration_stages_done:
                    st.pyplot(Flags(all_ms_files).plot_flags_occupancy())

            with st.expander("Direction Independent gains statistics", expanded=True):
                display_single_obsid_gains_stats(
                    self.nextleap_config_file, "mpi_di", "XX", "stations"
                )
                display_single_obsid_gains_stats(
                    self.nextleap_config_file, "mpi_di", "XX", "clusters"
                )

            with st.expander("Bandpass gains statistics", expanded=True):
                display_single_obsid_gains_stats(
                    self.nextleap_config_file, "bandpass", "XX", "stations"
                )
                display_single_obsid_gains_stats(
                    self.nextleap_config_file, "bandpass", "XX", "clusters"
                )

            with st.expander("Direction Dependent gains statistics", expanded=True):
                display_single_obsid_gains_stats(
                    self.nextleap_config_file, "mpi_dd", "XX", "stations"
                )
                display_single_obsid_gains_stats(
                    self.nextleap_config_file, "mpi_dd", "XX", "clusters"
                )

            with st.expander("All-sky images", expanded=True):
                self.plot_allsky_images()

            with st.expander("Power spectrum outputs", expanded=True):
                for pdf_name in self.pspipe_pdfs:
                    im = convert_pdf_to_jpeg(pdf_name)
                    st.image(im, caption=f"{str(pdf_name).strip('.pdf')}")

    def plot_allsky_images(self):
        self.di_fits_image = [
            self.di_fits_image[0]
            if self.di_fits_image and self._check_path_exists(self.di_fits_image[0])
            else None
        ][0]
        self.dd_fits_image = [
            self.dd_fits_image[0]
            if self.dd_fits_image and self._check_path_exists(self.dd_fits_image[0])
            else None
        ][0]

        # all sky maps display
        col1, col2 = st.columns(2)
        image_files = [self.di_fits_image, self.dd_fits_image]

        if all(image_files):
            figs = [all_sky_map(image_file) for image_file in image_files]
            col1.pyplot(figs[0])
            col2.pyplot(figs[1])
        else:
            if image_files[0] is not None:
                col1.caption(self.di_fits_image)
                col1.pyplot(all_sky_map(image_files[0]))
            else:
                col1.write("DI Image missing")
            if image_files[1] is not None:
                col2.caption(self.dd_fits_image)
                col2.pyplot(all_sky_map(image_files[1]))
            else:
                col2.write("DD Image missing")
        return

    def display_html_file(self, html_file):
        # hfile = open(html_file,'r')
        # content = hfile.read()
        # st.markdown(content,unsafe_allow_html=True)
        # hfile.close()
        p = open(html_file)
        components.html(p.read())
        return

    def display_trace_file(self, trace_file):
        trace = pd.read_csv(trace_file, sep="\t")
        st.dataframe(trace)
        return

    def _check_path_exists(self, file_path):
        return os.path.isfile(file_path)

    def all_ms_files(self, calibration_stage):
        # unravel the list of lists containing the subbands in each node
        all_ms_files = []
        for node in self.config["data"]["sub_bands_per_node"].keys():
            ms_files = [
                f"/net/{node}/{self.config['data']['path']}/{self.config[calibration_stage]['ms_pattern'].replace('???', str(sb_num).zfill(3))}"
                for sb_num in self.config["data"]["sub_bands_per_node"][node]
            ]
            all_ms_files += ms_files

        # subband_numbers = list(self.config["data"]["sub_bands_per_node"].values())
        # subband_numbers = [item for sublist in subband_numbers for item in sublist]

        # return the ms filenames involved in the requested calibration stage
        # ms_files = [
        #     {}{self.config[calibration_stage]["ms_pattern"].replace(
        #         "???", str(sb_num).zfill(3)
        #     )}
        #     for sb_num in subband_numbers
        # ]
        return all_ms_files


@st.cache_data
def build_gains(config_file, calibration_stage="mpi_di"):
    Obsid = ObsID(config_file)
    if calibration_stage in Obsid.calibration_stages_done:
        gs = Gains(Obsid, calibration=calibration_stage)
        gains_data = gs.get_all_gains()
        clusters = gs.get_cluster_ids()
        stations = gs.stations
        obsid = gs.obsid.obsid

        return obsid, stations, clusters, gains_data
    else:
        return -1


@st.cache_data
def display_single_obsid_gains_stats(config_file, calibration_stage, pol, stats_type):
    calibration_stage_gains_data = build_gains(
        config_file, calibration_stage=calibration_stage
    )
    if calibration_stage_gains_data == -1:
        st.write(f"{calibration_stage} calibration was not carried out!")
    else:
        obsid, stations, clusters, gains_data = calibration_stage_gains_data

        col1, col2, col3 = st.columns([1, 1, 1])
        if stats_type == "stations":
            med, rms_dt, rms_dnu = get_all_stations_amplitude_stats(
                gains_data, clusters, stations, pol=pol
            )
            with col1:
                fig = plot_gains_amplitude_stat_for_all_stations(
                    med,
                    "Gains median ampitude",
                    stations,
                    yrange=[np.log10(0.1), np.log10(10)],
                )
                st.plotly_chart(fig, use_container_width=True)
            with col2:
                fig = plot_gains_amplitude_stat_for_all_stations(
                    rms_dnu,
                    "Frequency RMS",
                    stations,
                    yrange=[np.log10(0.0001), np.log10(0.2)],
                )
                st.plotly_chart(fig, use_container_width=True)
            with col3:
                fig = plot_gains_amplitude_stat_for_all_stations(
                    rms_dt,
                    "Time RMS",
                    stations,
                    yrange=[np.log10(0.01), np.log10(10)],
                )
                st.plotly_chart(fig, use_container_width=True)

        if stats_type == "clusters":
            med, rms_dt, rms_dnu = get_all_clusters_amplitude_stats(
                gains_data, clusters, stations, pol=pol
            )

            with col1:
                fig = plot_gains_amplitude_stat_for_all_clusters(
                    med,
                    "Gains median ampitude",
                    clusters,
                    yrange=[np.log10(0.1), np.log10(10)],
                )
                st.plotly_chart(fig, use_container_width=True)
            with col2:
                fig = plot_gains_amplitude_stat_for_all_clusters(
                    rms_dnu,
                    "Frequency RMS",
                    clusters,
                    yrange=[np.log10(0.0001), np.log10(0.2)],
                )
                st.plotly_chart(fig, use_container_width=True)
            with col3:
                fig = plot_gains_amplitude_stat_for_all_clusters(
                    rms_dt,
                    "Time RMS",
                    clusters,
                    yrange=[np.log10(0.0001), np.log10(5)],
                )
                st.plotly_chart(fig, use_container_width=True)


@st.cache_data
def display_multiple_obsids_gains_stats(
    all_processed_nights: list[str],
    calibration_stage: str,
    polarization: str,
    statistic: str,  # either "stations" or  "clusters"
):
    for nightID in all_processed_nights:
        # with st.expander(f"{nightID}", expanded=True):
        display_single_obsid_gains_stats(
            Processed.nights[nightID], calibration_stage, polarization, statistic
        )
        st.caption(nightID)
