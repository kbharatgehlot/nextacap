import numpy as np
import streamlit as st
import plotly.graph_objects as go
from plotly.subplots import make_subplots


class Gains:
    def __init__(self, obsid: object, calibration="mpi_di"):
        self.obsid = obsid

        if calibration == "mpi_di":
            self.solutions_file = obsid.di_solutions_file
            self.eff_file = obsid.di_effective_clusters
            self.clusters_file = obsid.di_clusters_file

        if calibration == "bandpass":
            self.solutions_file = obsid.bandpass_solutions_file
            self.eff_file = obsid.bandpass_effective_clusters
            self.clusters_file = obsid.bandpass_clusters_file

        if calibration == "mpi_dd":
            self.solutions_file = obsid.dd_solutions_file
            self.eff_file = obsid.dd_effective_clusters
            self.clusters_file = obsid.dd_clusters_file


        self.meta_data = np.load(
            self.solutions_file.replace(".npy", ".npz").replace("_zsol", "")
        )

        self.stations = self.meta_data["stations"]

        self.freqs = self.meta_data["freqs"] * 1e-6

        self.eff_nr = np.load(self.eff_file)

        self.pols = dict(
            zip(["XX", "YY", "XY", "YX"], [[0, 0], [1, 1], [0, 1], [1, 0]])
        )
        self.pol_stokes = dict(
            zip(["I", "V", "U", "Q"], [[0, 0], [1, 1], [0, 1], [1, 0]])
        )

    def get_freqs(self):
        fmin, fmax = int(np.floor(self.freqs.min())), int(np.ceil(self.freqs.max()))
        return fmin, fmax

    def read_sagecal_solfile(self):
        data = np.load(self.solutions_file)
        data = data[:, :, :, :, :, 0] + 1j * data[:, :, :, :, :, 1]
        data = data.transpose(0, 1, 2, 4, 3)
        data = data.reshape(*(list(data.shape[:-1]) + [2, 2]))

        return data

    def get_cluster_ids(self):
        clusters = []
        with open(self.clusters_file) as f:
            for line in f.readlines():
                if line.strip().startswith("#"):
                    continue
                # print(line)
                s = line.split(" ")
                if s == ["\n"]:
                    continue
                clusters.append({"ID": int(s[0]), "NS": int(s[1]), "S": s[2:]})

        return dict(zip([v["ID"] for v in clusters], np.arange(len(clusters))))

    def get_gains(self, d, cluster, station, eff_nr):
        if cluster == 0:
            c_id = slice(None, int(eff_nr.cumsum()[cluster]))
        else:
            c_id = slice(
                int(eff_nr.cumsum()[cluster - 1]), int(eff_nr.cumsum()[cluster])
            )

        s = d[:, station, :, c_id]
        s = s.transpose(0, 2, 1, 3, 4)
        s = np.concatenate(s)

        return s

    def g_mul(self, g1, g2, c):
        return np.matmul(np.matmul(g1, c), g2.conj().transpose(0, 1, 3, 2))

    def cov2stokes(self, R):
        stokesI = 0.5 * (R[:, :, 0, 0] + R[:, :, 1, 1])
        stokesV = 0.5 * (-1j * (R[:, :, 0, 1] - R[:, :, 1, 0]))
        stokesQ = 0.5 * (R[:, :, 0, 0] - R[:, :, 1, 1])
        stokesU = 0.5 * (R[:, :, 0, 1] + R[:, :, 1, 0])

        return stokesI, stokesQ, stokesU, stokesV

    def get_all_gains(_self):
        """
        read sagecal calibration gains data from the output solutions file

        Parameters
        ----------
        stations : list
            The stations names
        clusters : range object/list
            Clusters ids
        d : array
            gains data
        eff_nr : _type_
            _description_

        Returns
        -------
        Dict
            A nested dictionary with the gains data per station per cluster.
            e.g.
            dict.keys: dict_keys([('CS001HBA0', 0), ('CS001HBA0', 1), ('CS001HBA0', 2), ('CS001HBA0', 3)...
            dict[('CS001HBA0', 0)].keys():  dict_keys(['XX', 'gg_XX', 'YY',
            'gg_YY', 'XY', 'gg_XY', 'YX', 'gg_YX', 'gg_I', 'gg_V', 'gg_U', 'gg_Q'])
        """
        d = _self.read_sagecal_solfile()

        gains = dict()
        for s_idx, s in enumerate(_self.stations):
            clusters = _self.get_cluster_ids()
            for c in range(len(clusters)):
                gains[(s, c)] = dict()
                # ref_station = 0
                # ref_g = _self.get_gains(d, c, ref_station, _self.eff_nr)
                g = _self.get_gains(d, c, s_idx, _self.eff_nr)  # / ref_g
                gg = _self.g_mul(g, g, np.eye(2))
                I, Q, U, V = _self.cov2stokes(gg)
                gg_stokes = np.stack(
                    [np.stack([I, Q], axis=-1), np.stack([U, V], axis=-1)], axis=-1
                )

                for pol_name, (i, j) in _self.pols.items():
                    gains[(s, c)][pol_name] = g[:, :, i, j]
                    gains[(s, c)]["gg_" + pol_name] = gg[:, :, i, j]
                for pol_name, (i, j) in _self.pol_stokes.items():
                    gains[(s, c)]["gg_" + pol_name] = gg_stokes[:, :, i, j]

        return gains


################################################################################################
# Gains plotting
###############################################################################################


@st.cache_data
def get_all_stations_amplitude_stats(gains_data, clusters, stations, pol="XX"):
    rms_dt = []
    med = []
    rms_dnu = []

    c_ii = range(len(clusters))

    for s_id, stn in enumerate(stations):
        rms_dt.append(
            [np.std(np.diff(gains_data[(stn, k)][pol], axis=0)) for k in c_ii]
        )
        med.append([np.median(abs(gains_data[(stn, k)][pol])) for k in c_ii])
        rms_dnu.append(
            [np.std(np.diff(gains_data[(stn, k)][pol], axis=1)) for k in c_ii]
        )

    rms_dt = np.array(rms_dt)
    med = np.array(med)
    rms_dnu = np.array(rms_dnu)

    return med, rms_dt, rms_dnu


@st.cache_data
def get_all_clusters_amplitude_stats(gains_data, clusters, stations, pol="XX"):
    rms_dt = []
    med = []
    rms_dnu = []

    for c in np.arange(0, len(clusters)):
        rms_dt.append(
            [np.std(np.diff(gains_data[(k, c)][pol], axis=0)) for k in stations]
        )
        med.append([np.median(abs(gains_data[(k, c)][pol])) for k in stations])
        rms_dnu.append(
            [np.std(np.diff(gains_data[(k, c)][pol], axis=1)) for k in stations]
        )

    rms_dt = np.array(rms_dt)
    med = np.array(med)
    rms_dnu = np.array(rms_dnu)

    return med, rms_dt, rms_dnu


@st.cache_data
def plot_gains_amplitude_stat_for_all_stations(stat, statname, stations, yrange):
    fig = make_subplots(
        rows=1,
        cols=1,
        # shared_xaxes=True,
        horizontal_spacing=1e-3,
        # vertical_spacing=1e-2,
    )  # , subplot_titles=[x[0] for x in list(pol_stokes.items())])
    for i, (stat, stat_name) in enumerate(zip([stat], [statname])):
        for s_id, stn in enumerate(stations):
            fig.add_trace(
                go.Box(
                    y=stat[s_id, :].T,
                    boxpoints=False,
                    name=f"{stn}",
                    jitter=0.3,
                    showlegend=False,
                    marker_color="rgb(7,40,89)",
                    line_color="rgb(7,40,89)",
                ),
                row=i + 1,
                col=1,
            )
            fig.update_yaxes(
                title=f"{stat_name}",
                linewidth=2,
                linecolor="gray",
                mirror=True,
                showgrid=True,
                griddash="dash",
                minor_griddash="dot",
                zeroline=False,
                row=i + 1,
                col=1,
                type="log",
                range=yrange,
            )
            fig.update_xaxes(
                linewidth=1,
                linecolor="black",
                mirror=True,
                showgrid=True,
                zeroline=False,
                tickvals=list(stations),
                ticktext=stations,
                row=i + 1,
                col=1,
            )
    fig.update_layout(
        template="plotly_white",
        margin=dict(l=20, r=20, t=50, b=20),  # title=pol
    )

    return fig


@st.cache_data
def plot_gains_amplitude_stat_for_all_clusters(stat, statname, clusters, yrange):
    c_ii = range(len(clusters))
    fig = make_subplots(
        rows=1,
        cols=1,
        # shared_xaxes=True,
        horizontal_spacing=1e-3,
        # vertical_spacing=1e-2,
    )  # , subplot_titles=[x[0] for x in list(pol_stokes.items())])
    for i, (stat, stat_name) in enumerate(zip([stat], [statname])):
        for cl in c_ii:
            fig.add_trace(
                go.Box(
                    y=stat[cl, :].T,
                    boxpoints=False,
                    name=f"clstr_{cl}",
                    jitter=0.3,
                    showlegend=False,
                    marker_color="rgb(7,40,89)",
                    line_color="rgb(7,40,89)",
                ),
                row=i + 1,
                col=1,
            )
            fig.update_yaxes(
                title=f"{stat_name}",
                linewidth=2,
                linecolor="gray",
                mirror=True,
                showgrid=True,
                griddash="dash",
                minor_griddash="dot",
                zeroline=False,
                row=i + 1,
                col=1,
                type="log",
                range=yrange,
            )
            fig.update_xaxes(
                linewidth=1,
                linecolor="black",
                mirror=True,
                showgrid=True,
                zeroline=False,
                tickvals=list(c_ii),
                ticktext=[f"cluster {c}" for c in list(clusters.keys())],
                row=i + 1,
                col=1,
            )

    fig.update_layout(
        template="plotly_white",
        margin=dict(l=20, r=20, t=50, b=20),  # title=pol
    )

    return fig
