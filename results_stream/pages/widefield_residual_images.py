import streamlit as st
from models.observations import ObsID, Processed

if st.button("refresh"):
    # Clear values from *all* all in-memory and on-disk data caches:
    st.cache_data.clear()


@st.cache_data
def display_all_obsids_skymaps(all_processed_observations: list[str]):
    with st.container():
        for observation in all_processed_observations:
            with st.expander(f"{observation}", expanded=True):
                Obsid = ObsID(Processed.nights[observation])
                Obsid.plot_allsky_images()


display_all_obsids_skymaps(list(Processed.nights.keys()))
