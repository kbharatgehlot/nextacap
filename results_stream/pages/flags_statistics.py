import streamlit as st
from models.observations import ObsID, Processed
from models.flags import Flags


if st.button("refresh flags"):
    # Clear values from *all* all in-memory and on-disk data caches:
    st.cache_data.clear()

cal = "mpi_dd"


@st.cache_data
def display_all_obsids_flags_occupancy(all_processed_observations: list[str]):
    with st.container():
        for observation in all_processed_observations:
            with st.expander(f"{observation}", expanded=True):
                Obsid = ObsID(Processed.nights[observation])
                filenames = Obsid.all_ms_files(cal)
                fig = Flags(filenames).plot_flags_occupancy(num_cores=20)
                st.pyplot(fig)


display_all_obsids_flags_occupancy(list(Processed.nights.keys()))
