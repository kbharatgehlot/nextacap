import numpy as np
import streamlit as st
import json
import os
from pathlib import Path
from .utils import convert_pdf_to_jpeg, all_sky_map

from . import gains


# @st.cache_resource
def mock_jsons(size):
    js = [
        "/home/users/chege/theleap/testdir/L254871/logs/params.json",
        # "/home/users/chege/theleap/testdir/L254871/logs/params.json",
        # "/home/users/chege/theleap/testdir/L254871/logs/params.json",
        # "/home/users/chege/theleap/testdir/L254871_v2/logs/params.json",
        "/home/users/chege/NCP/redshift1/L412388/logs/params.json",
    ]
    o = {}

    for i in range(size):
        n = "L2554" + str(i)
        o[n] = js[np.random.randint(0, 2)]
    return o


class Processed:
    nights = mock_jsons(3)


class ObsID:
    def __init__(self, nextleap_config_file):
        self.nextleap_config_file = nextleap_config_file

        with open(self.nextleap_config_file, "r") as f:
            self.config = json.load(f)

        self.obsid = self.config["data"]["obsid"]
        self.working_nodes = self.config["cluster"]["nodes"].split(",")
        self.master_node = self.config["cluster"]["masternode"]
        self.data_directory = (
            f"/net/{self.config['cluster']['masternode']}/{self.config['data']['path']}"
        )

        self.pspipe_pdfs = sorted(
            [path for path in Path(f"{self.data_directory}/ps").rglob("*.pdf")]
        )

        self.direction_independent_clusters_file = self.config["mpi_di"][
            "clusters_file"
        ]
        self.direction_dependent_clusters_file = self.config["mpi_dd"]["clusters_file"]

        self.direction_independent_effective_clusters = f"{self.data_directory}/solutions_sagecal_mpi_di/eff_nr_{os.path.basename(self.direction_independent_clusters_file)}.npy"
        self.direction_dependent_effective_clusters = f"{self.data_directory}/solutions_sagecal_mpi_dd/eff_nr_{os.path.basename(self.direction_dependent_clusters_file)}.npy"

        self.direction_independent_solutions_file = (
            f"{self.data_directory}/solutions_sagecal_mpi_di/{self.obsid}.npy"
        )
        self.bandpass_solutions_file = (
            f"{self.data_directory}/solutions_sagecal_bandpass/{self.obsid}.npy"
        )
        self.direction_dependent_solutions_file = (
            f"{self.data_directory}/solutions_sagecal_mpi_dd/{self.obsid}.npy"
        )

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
            with st.expander("NextLEAP configuration"):
                st.json(self.config)
            with st.expander("Direction Independent gains statistics", expanded=True):
                display_single_obsid_gains_stats(
                    self.nextleap_config_file, "DI", "XX", "stations"
                )
                display_single_obsid_gains_stats(
                    self.nextleap_config_file, "DI", "XX", "clusters"
                )

            with st.expander("Bandpass gains statistics", expanded=True):
                display_single_obsid_gains_stats(
                    self.nextleap_config_file, "BP", "XX", "stations"
                )
                display_single_obsid_gains_stats(
                    self.nextleap_config_file, "BP", "XX", "clusters"
                )

            with st.expander("Direction Dependent gains statistics", expanded=True):
                display_single_obsid_gains_stats(
                    self.nextleap_config_file, "DD", "XX", "stations"
                )
                display_single_obsid_gains_stats(
                    self.nextleap_config_file, "DD", "XX", "clusters"
                )

            with st.expander("All-sky images", expanded=True):
                self.plot_allsky_images()

            with st.expander("Power spectrum outputs", expanded=True):
                for pdf_name in self.pspipe_pdfs:
                    im = convert_pdf_to_jpeg(pdf_name)
                    st.image(im, caption=f"{str(pdf_name).strip('.pdf')}")

    def plot_allsky_images(self):
        di_image = f"{self.data_directory}/all_sky_DI_images/all_sky_DI-image.fits"
        dd_image = di_image.replace("DI", "DD")

        # all sky maps display
        col1, col2 = st.columns(2)

        image_files = [
            x if self._check_path_exists(x) else None for x in [di_image, dd_image]
        ]

        if all(image_files):
            figs = [all_sky_map(image_file) for image_file in image_files]
            col1.pyplot(figs[0])
            col2.pyplot(figs[1])
        else:
            if image_files[0] is not None:
                col1.pyplot(all_sky_map(image_files[0]))
            else:
                col1.write("DI Image missing")
            if image_files[1] is not None:
                col2.pyplot(all_sky_map(image_files[1]))
            else:
                col2.write("DD Image missing")
        return

    def _check_path_exists(self, file_path):
        return os.path.isfile(file_path)


@st.cache_data
def build(config_file, gains_type="DI"):
    Obsid = ObsID(config_file)
    gs = gains.Gains(Obsid, calibration=gains_type)
    gains_data = gs.get_all_gains()
    clusters = gs.get_cluster_ids()
    stations = gs.stations
    obsid = gs.obsid.obsid

    return obsid, stations, clusters, gains_data


@st.cache_data
def display_single_obsid_gains_stats(config_file, gains_type, pol, stats_type):
    obsid, stations, clusters, gains_data = build(config_file, gains_type=gains_type)

    col1, col2, col3 = st.columns([1, 1, 1])
    if stats_type == "stations":
        med, rms_dt, rms_dnu = gains.get_all_stations_amplitude_stats(
            gains_data, clusters, stations, pol=pol
        )
        with col1:
            fig = gains.plot_gains_amplitude_stat_for_all_stations(
                med,
                "Gains median ampitude",
                stations,
                yrange=[np.log10(0.1), np.log10(10)],
            )
            st.plotly_chart(fig, use_container_width=True)
        with col2:
            fig = gains.plot_gains_amplitude_stat_for_all_stations(
                rms_dnu,
                "Frequency RMS",
                stations,
                yrange=[np.log10(0.0001), np.log10(0.2)],
            )
            st.plotly_chart(fig, use_container_width=True)
        with col3:
            fig = gains.plot_gains_amplitude_stat_for_all_stations(
                rms_dt,
                "Time RMS",
                stations,
                yrange=[np.log10(0.01), np.log10(10)],
            )
            st.plotly_chart(fig, use_container_width=True)

    if stats_type == "clusters":
        med, rms_dt, rms_dnu = gains.get_all_clusters_amplitude_stats(
            gains_data, clusters, stations, pol=pol
        )

        with col1:
            fig = gains.plot_gains_amplitude_stat_for_all_clusters(
                med,
                "Gains median ampitude",
                clusters,
                yrange=[np.log10(0.1), np.log10(10)],
            )
            st.plotly_chart(fig, use_container_width=True)
        with col2:
            fig = gains.plot_gains_amplitude_stat_for_all_clusters(
                rms_dnu,
                "Frequency RMS",
                clusters,
                yrange=[np.log10(0.0001), np.log10(0.2)],
            )
            st.plotly_chart(fig, use_container_width=True)
        with col3:
            fig = gains.plot_gains_amplitude_stat_for_all_clusters(
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
