import streamlit as st
from models.observations import (
    Processed,
    display_multiple_obsids_gains_stats,
)


all_nights = list(Processed.nights.keys())

col1, col2, __, clear_cache = st.columns([0.15, 0.15, 0.55, 0.15])
with st.container():
    if col1.button("Station Statistics"):
        display_multiple_obsids_gains_stats(all_nights, "BP", "gg_I", "stations")

    if col2.button("Cluster Statistics"):
        display_multiple_obsids_gains_stats(all_nights, "BP", "gg_I", "clusters")

    if clear_cache.button("Clear BP cache"):
        # Clear values from *all* all in-memory and on-disk data caches:
        st.cache_data.clear()
